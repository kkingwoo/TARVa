"""
rnaseqpipe align: two-pass STAR alignment for RNA variant calling

    `tarva rnaseqpipe align --metadata <meta.csv> --qc-output <qc_dir>
        --output <aligned_dir> --ref <reference.fa[.gz]> --star-index <dir>
        [--gtf <annotation.gtf>] [--paired FWD,REV] [--sjdb-overhang N]
        [--threads N] [--jobs N] [--limit-bam-sort-ram BYTES]
        [--min-unique-pct F] [--sample-index N] [--overwrite] [--work-dir DIR]
        [--hpc] [--array] [--array-concurrency N] [--conda-env NAME]
        [--modules M1,M2] [--directives-dict "{...}"]`

This is the SECOND of the pipeline's two STAR alignments. The first
(`stringtiepipe align`) is run separately and those BAMs are created 
specifically for Stringtie.

Input is the trimmed RNA tree under --qc-output; output is
<output>/<group>/<individualID>/<individualID>_twopass.bam, which
rnaseqpipe/preprocess.py reads via TWOPASS_BAM_SUFFIX below.

"""

import csv
import subprocess
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tarva.metadata import DEFAULT_PAIRED_TOKENS, load_metadata, parse_paired_tokens
from tarva.reference import (DEFAULT_STAR_MODULE, ensure_star_index, ensure_unzipped,
                             ensure_unzipped_fasta)
from tarva.slurm_utils import parse_directives_dict, split_directives, submit_job, write_sbatch_script
from tarva.stringtiepipe.align import parse_star_log, rna_sjdb_overhang

TWOPASS_BAM_SUFFIX = "_twopass.bam"
DEFAULT_LIMIT_BAM_SORT_RAM = 30_000_000_000
DEFAULT_MIN_UNIQUE_PCT = 75.0
# Rough per-process genome-index footprint for a mammalian STAR index, used only
# for the --jobs memory warning.
MAMMALIAN_INDEX_GB = 30

SUMMARY_FIELDS = [
    "individualID", "group", "status", "bam", "input_reads",
    "uniquely_mapped_pct", "multi_mapped_pct", "too_short_pct",
    "splices_total", "error",
    ]


def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--qc-output", required=True,
                        help="The --output directory `tarva qc` wrote: supplies both the "
                             "trimmed/RNA/ tree and read_lengths.csv")
    parser.add_argument("--output", required=True,
                        help="Directory for two-pass BAMs, as <output>/<group>/<individualID>/")
    parser.add_argument("--ref", required=True,
                        help="Reference FASTA, decompressed first if it ends in .gz. Suggest using "
                             "combined_reference.fa from `stringtiepipe novel-transcripts`. ")
    parser.add_argument("--star-index", required=True, metavar="DIR",
                        help="STAR genome index directory. Must be the index built for --ref. "
                             "Built here if missing. ")
    parser.add_argument("--gtf",
                        help="Annotation GTF inserted at mapping time via --sjdbGTFfile. "
                             "OPTIONAL, unlike `stringtiepipe align`,  --twopassMode Basic "
                             "discovers junctions per sample on its own. See the "
                             "note in run() about combining this with a combined reference")
    parser.add_argument("--paired", type=parse_paired_tokens, default=DEFAULT_PAIRED_TOKENS,
                        metavar="FWD,REV",
                        help="Substrings identifying the forward and reverse mate in trimmed "
                             "FASTQ filenames (default: _1,_2), matched against the part of "
                             "the name following the metadata RNA prefix")
    parser.add_argument("--sjdb-overhang", type=int,
                        help="STAR --sjdbOverhang, only used with --gtf (default: derived "
                             "from read_lengths.csv as max RNA read length minus 1)")
    parser.add_argument("--work-dir",
                        help="Where to write the decompressed FASTA if --ref is gzipped "
                             "(default: alongside --ref)")
    parser.add_argument("--threads", type=int, default=8,
                        help="STAR --runThreadN per sample (default: 8)")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Samples to align concurrently within one process (default: 1). "
                             "Keep this at 1 unless you have the memory -- see "
                             "--limit-bam-sort-ram. Prefer --array on a cluster")
    parser.add_argument("--limit-bam-sort-ram", type=int, default=DEFAULT_LIMIT_BAM_SORT_RAM,
                        metavar="BYTES",
                        help=f"STAR --limitBAMsortRAM (default: {DEFAULT_LIMIT_BAM_SORT_RAM}), "
                             "Without this STAR fails sorting a large "
                             "coordinate-sorted BAM")
    parser.add_argument("--min-unique-pct", type=float, default=DEFAULT_MIN_UNIQUE_PCT,
                        help=f"Warn when a sample's uniquely-mapped %% falls below this "
                             f"(default: {DEFAULT_MIN_UNIQUE_PCT}). Samples are reported, "
                             "never dropped. Up to user discretion to manually remove those below "
                             "the specified threshold. ")
    parser.add_argument("--sample-index", type=int, metavar="N",
                        help="Align only sample N (0-based, metadata row order). Set "
                             "automatically from $SLURM_ARRAY_TASK_ID by --array")
    parser.add_argument("--overwrite", action="store_true",
                        help="Realign samples that already have a two-pass BAM")
    parser.add_argument("--hpc", action="store_true",
                        help="Submit as a SLURM job instead of running locally")
    parser.add_argument("--array", action="store_true",
                        help="With --hpc, submit a job array of one task per sample. Each "
                             "task gets the full --threads allocation and its own wall clock")
    parser.add_argument("--array-concurrency", type=int, metavar="N",
                        help="Cap simultaneously running array tasks (SLURM's %%N suffix). "
                             "Worth setting here: --genomeLoad NoSharedMemory means every "
                             "concurrent task holds its own copy of the genome index")
    parser.add_argument("--conda-env", default="bio_qc",
                        help="Conda environment activated in the generated SLURM script "
                             "(default: bio_qc -- runs the CLI itself and provides samtools)")
    parser.add_argument("--modules", default=DEFAULT_STAR_MODULE,
                        help=f"Comma-separated HPC modules to load (default: {DEFAULT_STAR_MODULE})")
    parser.add_argument("--directives-dict", dest="directives_dict",
                        type=parse_directives_dict, default={}, metavar="DICT",
                        help="SLURM sbatch directives as a Python dict literal, only used with "
                             "--hpc. Two-pass needs more wall clock than one pass. Example: "
                             "\"{'partition': 'Mamba', 'mem': '80gb', 'time': '24:00:00', "
                             "'cpus_per_task': 12}\"")


def _summary_row(sample, bam, sample_dir, status, error=""):
    stats = parse_star_log(Path(sample_dir) / "Log.final.out")
    return {
        "individualID": sample.individual_id,
        "group": sample.group,
        "status": status,
        "bam": str(bam) if bam else "",
        "input_reads": stats.get("Number of input reads", ""),
        "uniquely_mapped_pct": stats.get("Uniquely mapped reads %", "").rstrip("%"),
        "multi_mapped_pct": stats.get("% of reads mapped to multiple loci", "").rstrip("%"),
        "too_short_pct": stats.get("% of reads unmapped: too short", "").rstrip("%"),
        # Two-pass should raise this against the first pass on the same sample
        "splices_total": stats.get("Number of splices: Total", ""),
        "error": error,
        }


def _align_one(sample, qc_output, output_dir, index_dir, gtf, sjdb_overhang,
               paired_tokens, threads, limit_bam_sort_ram, overwrite):
    """
    One sample's two-pass alignment. Failures become FAILED rows rather than
    exceptions, so one bad sample doesn't abort a pool with hours in flight.
    """
    sample_dir = Path(output_dir) / sample.group / sample.individual_id
    final_bam = sample_dir / f"{sample.individual_id}{TWOPASS_BAM_SUFFIX}"
    try:
        if final_bam.is_file() and not overwrite:
            print(f"{sample.individual_id}: {final_bam.name} already present, skipping "
                  f"(pass --overwrite to realign)")
            return _summary_row(sample, final_bam, sample_dir, status="skipped")

        fwd, rev = sample.trimmed_rna_pair(qc_output, paired_tokens)
        sample_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            "STAR",
            "--runThreadN", str(threads),
            "--genomeDir", str(index_dir),
            # Mandatory with --twopassMode, and with on-the-fly --sjdbGTFfile
            # insertion. Stale shared segments from a LoadAndKeep run are a
            # recurring source of unrelated-looking failures.
            "--genomeLoad", "NoSharedMemory",
            "--twopassMode", "Basic",
            "--readFilesCommand", "zcat",
            "--readFilesIn", str(fwd), str(rev),
            "--outFileNamePrefix", f"{sample_dir}/",
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outSAMmode", "Full",
            # NM and MD are the mismatch tags the GATK tools expect.
            "--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD",
            "--limitBAMsortRAM", str(limit_bam_sort_ram),
            ]
        if gtf:
            cmd += ["--sjdbGTFfile", str(gtf), "--sjdbOverhang", str(sjdb_overhang)]
        print(f"{sample.individual_id}: " + " ".join(cmd), flush=True)
        result = subprocess.run(cmd)
        if result.returncode != 0:
            return _summary_row(
                    sample, None, sample_dir, status="FAILED",
                    error=f"STAR exited {result.returncode}; see {sample_dir / 'Log.out'}",
                    )

        star_bam = sample_dir / "Aligned.sortedByCoord.out.bam"
        if not star_bam.is_file():
            return _summary_row(
                    sample, None, sample_dir, status="FAILED",
                    error=f"STAR exited 0 but {star_bam.name} is missing",
                    )
        star_bam.replace(final_bam)
        return _summary_row(sample, final_bam, sample_dir, status="aligned")
    except Exception as exc:
        return _summary_row(sample, None, sample_dir, status="FAILED",
                            error=f"{type(exc).__name__}: {exc}")


def _select_samples(meta, sample_index):
    if sample_index is None:
        return meta.samples
    if not 0 <= sample_index < len(meta.samples):
        raise SystemExit(
                f"rnaseqpipe align: --sample-index {sample_index} is out of range for the "
                f"{len(meta.samples)} sample(s) in --metadata "
                f"(valid: 0..{len(meta.samples) - 1})."
                )
    return [meta.samples[sample_index]]


def _summary_path(output_dir, samples, sample_index):
    """."""
    if sample_index is None:
        return Path(output_dir) / "twopass_alignment_summary.csv"
    return (Path(output_dir) / "summary_parts"
            / f"{sample_index:03d}_{samples[0].individual_id}.csv")


def _merge_summary_parts(output_dir, total_samples):
    """Last array task to finish assembles the summary from the parts."""
    parts = sorted((Path(output_dir) / "summary_parts").glob("*.csv"))
    combined = Path(output_dir) / "twopass_alignment_summary.csv"
    if len(parts) < total_samples:
        print(f"Summary parts: {len(parts)}/{total_samples}. {combined.name} will be "
              f"written by the last array task to finish.")
        return None
    partial = combined.with_name(combined.name + ".partial")
    with open(partial, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for part in parts:
            with open(part, newline="") as part_fh:
                for row in csv.DictReader(part_fh):
                    writer.writerow(row)
    partial.replace(combined)
    print(f"Wrote combined summary: {combined}")
    return combined


def run_alignment(args, samples, index_dir, gtf, sjdb_overhang, total_samples):
    """Across --jobs samples, flush the summary per completed sample."""
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = _summary_path(output_dir, samples, args.sample_index)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    rows = []

    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        fh.flush()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [
                pool.submit(_align_one, sample, args.qc_output, output_dir, index_dir,
                            gtf, sjdb_overhang, args.paired, args.threads,
                            args.limit_bam_sort_ram, args.overwrite)
                for sample in samples
                ]
            for future in as_completed(futures):
                row = future.result()
                writer.writerow(row)
                fh.flush()
                rows.append(row)

    aligned = [r for r in rows if r["status"] == "aligned"]
    skipped = [r for r in rows if r["status"] == "skipped"]
    failed = [r for r in rows if r["status"] == "FAILED"]
    print(f"\nTwo-pass alignment: {len(aligned)} aligned, {len(skipped)} skipped, "
          f"{len(failed)} failed")
    print(f"Wrote summary: {summary_csv}")
    if args.sample_index is not None:
        _merge_summary_parts(output_dir, total_samples)

    for row in failed:
        print(f"  FAILED {row['individualID']}: {row['error']}", file=sys.stderr)
    for row in rows:
        pct = row["uniquely_mapped_pct"]
        if pct and float(pct) < args.min_unique_pct:
            print(f"WARNING: {row['individualID']} uniquely mapped {pct}%, below "
                  f"--min-unique-pct {args.min_unique_pct}. Reports reads, does not remove them "
                  f"User descretion to exclude reads deliberately.",
                  file=sys.stderr)
    if failed:
        raise SystemExit(
                f"rnaseqpipe align: {len(failed)} sample(s) failed -- see {summary_csv}"
                )
    return rows


def _warn_memory(args):
    """
    NoSharedMemory means each concurrent STAR holds its own genome copy, and
    --twopassMode adds the re-inserted junction database on top.
    """
    if args.jobs <= 1:
        return
    sort_gb = args.limit_bam_sort_ram / 1e9
    print(f"WARNING: --jobs {args.jobs} runs {args.jobs} STAR processes at once, and "
          f"--genomeLoad NoSharedMemory gives each its own copy of the genome index "
          f"(~{MAMMALIAN_INDEX_GB}gb for a mammalian genome) plus up to {sort_gb:.0f}gb of "
          f"BAM sort buffer, before --twopassMode's re-inserted junctions. Budget at "
          f"least {args.jobs * (MAMMALIAN_INDEX_GB + sort_gb):.0f}gb of --mem, or use "
          f"--jobs 1 with --array.", file=sys.stderr)


def run(args):
    """
    Note on --gtf with a combined reference:

    combined_reference.fa carries the novel transcripts as their own contigs,
    while merged.gtf annotates those same transcripts as features on the genomic
    contigs.
    """
    if args.hpc:
        _submit_hpc(args)
        return
    meta = load_metadata(args.metadata)
    samples = _select_samples(meta, args.sample_index)
    fasta = ensure_unzipped_fasta(args.ref, work_dir=args.work_dir)
    gtf = (ensure_unzipped(args.gtf, work_dir=args.work_dir, description="annotation GTF")
           if args.gtf else None)
    sjdb_overhang = None
    if gtf:
        sjdb_overhang = args.sjdb_overhang or rna_sjdb_overhang(args.qc_output)
    index_dir = ensure_star_index(fasta, args.star_index, threads=args.threads)

    print(f"Reference: {fasta}")
    print(f"Two-pass STAR: {len(samples)} of {len(meta.samples)} sample(s) "
          f"({', '.join(s.individual_id for s in samples)}), {args.jobs} at a time, "
          f"{args.threads} threads each")
    print(f"Annotation: {gtf if gtf else 'none (--twopassMode Basic discovers junctions)'}"
          + (f", --sjdbOverhang {sjdb_overhang}" if gtf else ""))
    _warn_memory(args)
    run_alignment(args, samples, index_dir, gtf, sjdb_overhang, len(meta.samples))


def _submit_hpc(args):
    command_parts = [
        # -u: unbuffered, so a job killed at the wall clock still leaves its
        # progress in the .o rather than discarding the buffer.
        "python", "-u", "-m", "tarva.cli", "rnaseqpipe", "align",
        "--metadata", str(args.metadata),
        "--qc-output", str(args.qc_output),
        "--output", str(args.output),
        "--ref", str(args.ref),
        "--star-index", str(args.star_index),
        "--paired", ",".join(args.paired),
        "--threads", str(args.threads),
        "--jobs", str(args.jobs),
        "--limit-bam-sort-ram", str(args.limit_bam_sort_ram),
        "--min-unique-pct", str(args.min_unique_pct),
        ]
    for flag, value in (("--gtf", args.gtf), ("--sjdb-overhang", args.sjdb_overhang),
                        ("--work-dir", args.work_dir)):
        if value is not None:
            command_parts += [flag, str(value)]
    if args.overwrite:
        command_parts.append("--overwrite")

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    sbatch_kwargs.setdefault("cpus_per_task", args.threads * args.jobs)

    log_suffix = ""
    if args.array:
        if args.jobs > 1:
            print(f"WARNING: --array already gives each sample its own task, so --jobs "
                  f"{args.jobs} multiplies every task's genome-index memory. Use "
                  f"--jobs 1 with --array.", file=sys.stderr)
        n_samples = len(load_metadata(args.metadata).samples)
        array_spec = f"0-{n_samples - 1}"
        if args.array_concurrency:
            array_spec += f"%{args.array_concurrency}"
        extra_directives = dict(extra_directives)
        extra_directives.setdefault("array", array_spec)
        command_parts += ["--sample-index", "$SLURM_ARRAY_TASK_ID"]
        log_suffix = "_%a"
        print(f"Submitting a job array of {n_samples} task(s): --array={array_spec}")

    log_dir = Path(args.output) / "logs"
    script_path = write_sbatch_script(
            job_name="tarva_rna_align",
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / "tarva_rna_align.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            log_suffix=log_suffix,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
