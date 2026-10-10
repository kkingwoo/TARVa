"""
stringtiepipe align: first-pass STAR alignment, StringTie's prerequisite

    `tarva stringtiepipe align --metadata <meta.csv> --qc-output <qc_dir>
        --output <aligned_dir> --ref <genome.fa[.gz]> --star-index <dir>
        --gtf <annotation.gtf> [--paired FWD,REV] [--sjdb-overhang N]
        [--threads N] [--jobs N] [--limit-bam-sort-ram BYTES]
        [--min-unique-pct F] [--overwrite] [--work-dir DIR]
        [--hpc] [--conda-env NAME] [--modules M1,M2] [--directives-dict "{...}"]`

This is plain STAR, deliberately WITHOUT --twopassMode: StringTie needs an
aligned BAM and cannot run on FASTQ, so this pass exists only to feed it.
The separate two-pass alignment for GATK variant calling is
rnaseqpipe/align.py and is run later, after StringTie finishes.

Input is the trimmed RNA tree under --qc-output (see metadata.py); output is
<output>/<group>/<individualID>/, with STAR's Aligned.sortedByCoord.out.bam
renamed to <individualID>_firstpass.bam so assemble.py has a predictable path.
STAR's own logs (Log.final.out, SJ.out.tab) are left in place beside it.

A per-sample failure is logged and the remaining samples continue, but the
command exits non-zero if any sample failed, so a chained --hpc job using
'dependency': 'afterok:<id>' won't proceed on partial results.
"""

import csv
import subprocess
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from tarva.metadata import DEFAULT_PAIRED_TOKENS, load_metadata, parse_paired_tokens
from tarva.reference import (DEFAULT_STAR_MODULE, ensure_star_index, ensure_unzipped,
                             ensure_unzipped_fasta)
from tarva.slurm_utils import parse_directives_dict, split_directives, submit_job, write_sbatch_script

DEFAULT_LIMIT_BAM_SORT_RAM = 30_000_000_000
DEFAULT_MIN_UNIQUE_PCT = 75.0
FIRSTPASS_BAM_SUFFIX = "_firstpass.bam"
# Rough per-process genome-index footprint for a mammalian STAR index, used
# only for the --jobs memory warning below.
MAMMALIAN_INDEX_GB = 30

SUMMARY_FIELDS = [
    "individualID", "group", "status", "bam", "input_reads",
    "uniquely_mapped_pct", "multi_mapped_pct", "too_short_pct", "error",
    ]


def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--qc-output", required=True,
                        help="The --output directory `tarva qc` wrote: supplies both the "
                             "trimmed/ tree and read_lengths.csv")
    parser.add_argument("--output", required=True,
                        help="Directory for first-pass BAMs, as <output>/<group>/<individualID>/")
    parser.add_argument("--ref", required=True,
                        help="Reference genome FASTA; decompressed first if it ends in .gz")
    parser.add_argument("--star-index", required=True, metavar="DIR",
                        help="STAR genome index directory (built by `tarva prep-reference`; "
                             "built here if missing)")
    parser.add_argument("--gtf", required=True,
                        help="Annotation GTF, inserted at mapping time via --sjdbGTFfile")
    parser.add_argument("--paired", type=parse_paired_tokens, default=DEFAULT_PAIRED_TOKENS,
                        metavar="FWD,REV",
                        help="Substrings identifying the forward and reverse mate in trimmed "
                             "FASTQ filenames (default: _1,_2), matched against the part of "
                             "the name following the metadata RNA prefix")
    parser.add_argument("--sjdb-overhang", type=int,
                        help="STAR --sjdbOverhang (default: derived from read_lengths.csv as "
                             "max RNA read length minus 1)")
    parser.add_argument("--work-dir",
                        help="Where to write the decompressed FASTA if --ref is gzipped "
                             "(default: alongside --ref)")
    parser.add_argument("--threads", type=int, default=8,
                        help="STAR --runThreadN per sample (default: 8)")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Samples to align concurrently (default: 1). Keep this at 1 "
                             "unless you have the memory for it -- see --limit-bam-sort-ram")
    parser.add_argument("--limit-bam-sort-ram", type=int, default=DEFAULT_LIMIT_BAM_SORT_RAM,
                        metavar="BYTES",
                        help=f"STAR --limitBAMsortRAM (default: {DEFAULT_LIMIT_BAM_SORT_RAM}). "
                             "Without this STAR fails sorting a large coordinate-sorted BAM. "
                             "Each concurrent job needs its own, on top of its own genome copy")
    parser.add_argument("--min-unique-pct", type=float, default=DEFAULT_MIN_UNIQUE_PCT,
                        help=f"Warn when a sample's uniquely-mapped %% falls below this "
                             f"(default: {DEFAULT_MIN_UNIQUE_PCT}). Samples are reported, "
                             "never dropped")
    parser.add_argument("--overwrite", action="store_true",
                        help="Realign samples that already have a first-pass BAM")
    parser.add_argument("--hpc", action="store_true",
                        help="Submit as a SLURM job instead of running locally")
    parser.add_argument("--conda-env", default="bio_qc",
                        help="Conda environment activated in the generated SLURM script "
                             "(default: bio_qc)")
    parser.add_argument("--modules", default=DEFAULT_STAR_MODULE,
                        help=f"Comma-separated HPC modules to load (default: {DEFAULT_STAR_MODULE})")
    parser.add_argument("--directives-dict", dest="directives_dict",
                        type=parse_directives_dict, default={}, metavar="DICT",
                        help="SLURM sbatch directives as a Python dict literal, only used with "
                             "--hpc. Example: \"{'partition': 'Mambaa', 'mem': '80gb', "
                             "'time': '24:00:00', 'cpus_per_task': 8}\"")


def rna_sjdb_overhang(qc_output):
    """
    STAR --sjdbOverhang, which the manual defines as max read length - 1.

    """
    csv_path = Path(qc_output) / "read_lengths.csv"
    if not csv_path.is_file():
        raise SystemExit(
                f"stringtiepipe align: {csv_path} not found, so --sjdb-overhang cannot be "
                f"derived. Point --qc-output at the directory `tarva qc` wrote, or pass "
                f"--sjdb-overhang explicitly (max RNA read length minus 1)."
                )
    df = pd.read_csv(csv_path)
    rna = df[df["seq_type"] == "RNA"]
    if rna.empty:
        raise SystemExit(f"stringtiepipe align: no RNA rows in {csv_path}")
    return int(rna["max_length"].max()) - 1


def parse_star_log(log_path):
    """Flatten STAR's Log.final.out 'name | value' table into a dict."""
    stats = {}
    log_path = Path(log_path)
    if not log_path.is_file():
        return stats
    with open(log_path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            key, _, value = line.partition("|")
            stats[key.strip()] = value.strip()
    return stats


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
        "error": error,
        }


def _align_one(sample, qc_output, output_dir, index_dir, gtf, sjdb_overhang,
               paired_tokens, threads, limit_bam_sort_ram, overwrite):
    """
    First-pass alignment for a single individual.

    """
    sample_dir = Path(output_dir) / sample.group / sample.individual_id
    final_bam = sample_dir / f"{sample.individual_id}{FIRSTPASS_BAM_SUFFIX}"
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
            # Mandatory here: on-the-fly --sjdbGTFfile insertion is incompatible
            # with shared-memory genome loading, and stale shared segments from a
            # LoadAndKeep run are a recurring source of unrelated-looking failures.
            "--genomeLoad", "NoSharedMemory",
            "--sjdbGTFfile", str(gtf),
            "--sjdbOverhang", str(sjdb_overhang),
            "--readFilesCommand", "zcat",
            "--readFilesIn", str(fwd), str(rev),
            "--outFileNamePrefix", f"{sample_dir}/",
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outSAMstrandField", "intronMotif",
            # StringTie requires IH/HI attributes to start at 0.
            "--outSAMattrIHstart", "0",
            "--chimOutType", "WithinBAM",
            "--limitBAMsortRAM", str(limit_bam_sort_ram),
            ]
        print(f"{sample.individual_id}: " + " ".join(cmd))
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
        return _summary_row(sample, None, sample_dir, status="FAILED", error=f"{type(exc).__name__}: {exc}")


def run_alignment(args, meta, index_dir, gtf, sjdb_overhang):
    """ 
    """
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = output_dir / "firstpass_alignment_summary.csv"
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
                for sample in meta.samples
                ]
            for future in as_completed(futures):
                row = future.result()
                writer.writerow(row)
                fh.flush()
                rows.append(row)

    aligned = [r for r in rows if r["status"] == "aligned"]
    skipped = [r for r in rows if r["status"] == "skipped"]
    failed = [r for r in rows if r["status"] == "FAILED"]
    print(f"\nFirst-pass alignment: {len(aligned)} aligned, {len(skipped)} skipped, "
          f"{len(failed)} failed")
    print(f"Wrote summary: {summary_csv}")

    for row in failed:
        print(f"  FAILED {row['individualID']}: {row['error']}", file=sys.stderr)
    for row in rows:
        pct = row["uniquely_mapped_pct"]
        if pct and float(pct) < args.min_unique_pct:
            print(f"WARNING: {row['individualID']} uniquely mapped {pct}%, below "
                  f"--min-unique-pct {args.min_unique_pct}. Reads falling under threshold "
                  f"are reported, not dropped. Remove them separately if needed.",
                  file=sys.stderr)
    if failed:
        raise SystemExit(
                f"stringtiepipe align: {len(failed)} sample(s) failed -- see {summary_csv}"
                )
    return rows


def _warn_memory(args):
    """
    NoSharedMemory means each concurrent STAR holds its own genome copy, so
    memory scales with --jobs.

    """
    if args.jobs <= 1:
        return
    sort_gb = args.limit_bam_sort_ram / 1e9
    print(f"WARNING: --jobs {args.jobs} runs {args.jobs} STAR processes at once, and "
          f"--genomeLoad NoSharedMemory gives each its own copy of the genome index "
          f"(~{MAMMALIAN_INDEX_GB}gb for a mammalian genome) plus up to {sort_gb:.0f}gb of "
          f"BAM sort buffer. Budget roughly "
          f"{args.jobs * (MAMMALIAN_INDEX_GB + sort_gb):.0f}gb of --mem, or use --jobs 1.",
          file=sys.stderr)


def run(args):
    if args.hpc:
        _submit_hpc(args)
        return
    meta = load_metadata(args.metadata)
    sjdb_overhang = args.sjdb_overhang or rna_sjdb_overhang(args.qc_output)
    fasta = ensure_unzipped_fasta(args.ref, work_dir=args.work_dir)
    gtf = ensure_unzipped(args.gtf, work_dir=args.work_dir, description="annotation GTF")
    index_dir = ensure_star_index(fasta, args.star_index, threads=args.threads)
    print(f"First-pass STAR: {len(meta.samples)} sample(s), --sjdbOverhang {sjdb_overhang}, "
          f"{args.jobs} at a time, {args.threads} threads each")
    _warn_memory(args)
    run_alignment(args, meta, index_dir, gtf, sjdb_overhang)


def _submit_hpc(args):
    command_parts = [
        # -u: unbuffered, so a job killed at the wall clock still leaves its
        # progress in the .o rather than discarding the buffer.
        "python", "-u", "-m", "tarva.cli", "stringtiepipe", "align",
        "--metadata", str(args.metadata),
        "--qc-output", str(args.qc_output),
        "--output", str(args.output),
        "--ref", str(args.ref),
        "--star-index", str(args.star_index),
        "--gtf", str(args.gtf),
        "--paired", ",".join(args.paired),
        "--threads", str(args.threads),
        "--jobs", str(args.jobs),
        "--limit-bam-sort-ram", str(args.limit_bam_sort_ram),
        "--min-unique-pct", str(args.min_unique_pct),
        ]
    if args.sjdb_overhang is not None:
        command_parts += ["--sjdb-overhang", str(args.sjdb_overhang)]
    if args.work_dir:
        command_parts += ["--work-dir", str(args.work_dir)]
    if args.overwrite:
        command_parts.append("--overwrite")

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    sbatch_kwargs.setdefault("cpus_per_task", args.threads * args.jobs)
    log_dir = Path(args.output) / "logs"
    script_path = write_sbatch_script(
            job_name="tarva_stringtie_align",
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / "tarva_stringtie_align.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
