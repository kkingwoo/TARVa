"""
wgseqpipe align: BWA-mem alignment of trimmed WGS reads

    `tarva wgseqpipe align --metadata <meta.csv> --qc-output <qc_dir>
        --output <aligned_dir> --ref <genome.fa[.gz]> [--paired FWD,REV]
        [--threads N] [--sort-threads N] [--sort-memory SIZE] [--temp-dir DIR]
        [--jobs N] [--overwrite] [--work-dir DIR] [--build-index]
        [--sample-index N]
        [--hpc] [--array] [--array-concurrency N] [--conda-env NAME]
        [--modules M1,M2] [--directives-dict "{...}"]`

BWA-mem pipes straight into samtools sort, with no intermediate SAM, and both
return codes are checked -- bwa dying mid-stream otherwise leaves a
truncated-but-valid BAM and exit 0.

Input:  the trimmed WGS tree under --qc-output (see metadata.py) and a
        reference with its BWA index already beside it (--build-index opts
        into building a missing one, which takes 1-2 hours single-threaded).
Output: <output>/<group>/<individualID>/<individualID>_wgs.bam,
        coordinate-sorted and indexed, plus a summary CSV flushed per sample.
        --array gives each sample its own SLURM task and its own wall clock.
"""

import csv
import os
import subprocess
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tarva.metadata import DEFAULT_PAIRED_TOKENS, load_metadata, parse_paired_tokens
from tarva.reference import DEFAULT_BWA_MODULE, ensure_bwa_index, ensure_unzipped_fasta
from tarva.slurm_utils import parse_directives_dict, split_directives, submit_job, write_sbatch_script

ALIGNED_BAM_SUFFIX = "_wgs.bam"
DEFAULT_SORT_MEMORY = "2G"
BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")

SUMMARY_FIELDS = [
    "individualID", "group", "status", "bam", "mapped_reads", "unmapped_reads", "error",
    ]


def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--qc-output", required=True,
                        help="The --output directory `tarva qc` wrote; trimmed WGS reads are "
                             "read from its trimmed/WGS/ tree")
    parser.add_argument("--output", required=True,
                        help="Directory for aligned BAMs, as <output>/<group>/<individualID>/")
    parser.add_argument("--ref", required=True,
                        help="Reference genome FASTA; decompressed first if it ends in .gz. "
                             "The BWA index must already sit beside it (see --build-index)")
    parser.add_argument("--paired", type=parse_paired_tokens, default=DEFAULT_PAIRED_TOKENS,
                        metavar="FWD,REV",
                        help="Substrings identifying the forward and reverse mate in trimmed "
                             "FASTQ filenames (default: _1,_2), matched against the part of "
                             "the name following the metadata WGS prefix")
    parser.add_argument("--work-dir",
                        help="Where to write the decompressed FASTA if --ref is gzipped "
                             "(default: alongside --ref)")
    parser.add_argument("--temp-dir", metavar="DIR",
                        help="Root for samtools sort temporary files (default: the sample's "
                             "own output directory). ")
    parser.add_argument("--build-index", action="store_true",
                        help="Build the BWA index here if it is missing. Off by default: "
                             "`tarva prep-reference --bwa-index` ")
    parser.add_argument("--threads", type=int, default=8,
                        help="bwa mem -t per sample (default: 8)")
    parser.add_argument("--sort-threads", type=int, default=4,
                        help="samtools sort -@ per sample (default: 4). These run concurrently "
                             "with bwa in the same pipe, so budget --threads + --sort-threads "
                             "cores per job")
    parser.add_argument("--sort-memory", default=DEFAULT_SORT_MEMORY, metavar="SIZE",
                        help=f"samtools sort -m, per sort thread (default: {DEFAULT_SORT_MEMORY}). "
                             "Total sort memory is roughly this times --sort-threads")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Samples to align concurrently within one process (default: 1). "
                             "Each holds its own copy of the BWA index (~5gb for a mammalian "
                             "genome) plus its sort buffers. Prefer --array on a cluster")
    parser.add_argument("--sample-index", type=int, metavar="N",
                        help="Align only sample N (0-based, metadata row order) instead of all "
                             "of them. Set automatically from $SLURM_ARRAY_TASK_ID by --array; "
                             "rarely useful by hand")
    parser.add_argument("--overwrite", action="store_true",
                        help="Realign samples that already have a BAM")
    parser.add_argument("--hpc", action="store_true",
                        help="Submit as a SLURM job instead of running locally")
    parser.add_argument("--array", action="store_true",
                        help="With --hpc, submit a SLURM job array of one task per sample ")
    parser.add_argument("--array-concurrency", type=int, metavar="N",
                        help="Cap simultaneously running array tasks (SLURM's %%N suffix). "
                             "Default is uncapped")
    parser.add_argument("--conda-env", default="bio_qc",
                        help="Conda environment activated in the generated SLURM script "
                             "(default: bio_qc -- runs the CLI itself and provides samtools)")
    parser.add_argument("--modules", default=DEFAULT_BWA_MODULE,
                        help=f"Comma-separated HPC modules to load (default: {DEFAULT_BWA_MODULE})")
    parser.add_argument("--directives-dict", dest="directives_dict",
                        type=parse_directives_dict, default={}, metavar="DICT",
                        help="SLURM sbatch directives as a Python dict literal, only used with "
                             "--hpc. WGS alignment is slow -- budget many hours per sample even "
                             "with --array. Example: \"{'partition': 'Mamba', 'mem': '32gb', "
                             "'time': '48:00:00', 'cpus_per_task': 16}\"")


def _sample_dir(root, sample):
    return Path(root) / sample.group / sample.individual_id


def _resolve_bwa_index(fasta, build_index):
    """
    Confirm bwa mem's index exists rather than quietly building it.

    """
    fasta = Path(fasta)
    missing = [s for s in BWA_INDEX_SUFFIXES
               if not fasta.with_name(fasta.name + s).is_file()]
    if not missing:
        print(f"BWA index found beside {fasta}")
        return fasta
    if not build_index:
        raise SystemExit(
                f"wgseqpipe align: no BWA index beside {fasta} -- missing {missing}. "
                f"Run `tarva prep-reference --ref {fasta} --bwa-index` ")
    print(f"--build-index: {missing} missing beside {fasta}, building now "
          f"(single-threaded)")
    return ensure_bwa_index(fasta)


def _select_samples(meta, sample_index):
    """All samples, or just the one an array task was handed."""
    if sample_index is None:
        return meta.samples
    if not 0 <= sample_index < len(meta.samples):
        raise SystemExit(
                f"wgseqpipe align: --sample-index {sample_index} is out of range for the "
                f"{len(meta.samples)} sample(s) in --metadata (valid: 0..{len(meta.samples) - 1}). "
                )
    return [meta.samples[sample_index]]


def _summary_path(output_dir, samples, sample_index):
    """
    One writer, one file.

    """
    if sample_index is None:
        return Path(output_dir) / "wgs_alignment_summary.csv"
    return (Path(output_dir) / "summary_parts"
            / f"{sample_index:03d}_{samples[0].individual_id}.csv")


def _merge_summary_parts(output_dir, total_samples):
    """
    Combine array parts into wgs_alignment_summary.csv once every part exists.

    """
    parts = sorted((Path(output_dir) / "summary_parts").glob("*.csv"))
    combined = Path(output_dir) / "wgs_alignment_summary.csv"
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
    os.replace(partial, combined)
    print(f"Wrote combined summary: {combined}")
    return combined


def _idxstats(bam):
    """
    Mapped/unmapped totals from the BAM index.

    """
    try:
        result = subprocess.run(["samtools", "idxstats", str(bam)],
                                capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, OSError):
        return "", ""
    mapped = unmapped = 0
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        mapped += int(fields[2])
        unmapped += int(fields[3])
    return mapped, unmapped


def _summary_row(sample, bam, status, error=""):
    mapped, unmapped = _idxstats(bam) if bam and Path(bam).is_file() else ("", "")
    return {
        "individualID": sample.individual_id,
        "group": sample.group,
        "status": status,
        "bam": str(bam) if bam else "",
        "mapped_reads": mapped,
        "unmapped_reads": unmapped,
        "error": error,
        }


def _bwa_sort(ref, fwd, rev, out_bam, temp_prefix, threads, sort_threads, sort_memory):
    """
    bwa mem piped into samtools sort, returning (bwa_rc, sort_rc).

    """
    bwa_cmd = ["bwa", "mem", "-t", str(threads), str(ref), str(fwd), str(rev)]
    sort_cmd = ["samtools", "sort", "-@", str(sort_threads), "-m", str(sort_memory),
                "-T", str(temp_prefix), "-o", str(out_bam), "-"]
    print("Running: " + " ".join(bwa_cmd) + " | " + " ".join(sort_cmd))
    bwa = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
    sort = subprocess.Popen(sort_cmd, stdin=bwa.stdout)
    bwa.stdout.close()
    sort_rc = sort.wait()
    bwa_rc = bwa.wait()
    return bwa_rc, sort_rc


def _clear_sort_temps(temp_prefix):
    """
    Remove <prefix>.NNNN.bam left behind by a killed sort.

    """
    temp_prefix = Path(temp_prefix)
    for leftover in temp_prefix.parent.glob(f"{temp_prefix.name}.*.bam"):
        try:
            leftover.unlink()
        except OSError:
            pass


def _align_one(sample, qc_output, output_dir, ref, temp_dir, paired_tokens, threads,
               sort_threads, sort_memory, overwrite):
    """
    One sample's alignment. 
    """
    sample_dir = _sample_dir(output_dir, sample)
    final_bam = sample_dir / f"{sample.individual_id}{ALIGNED_BAM_SUFFIX}"
    partial_bam = final_bam.with_name(final_bam.name + ".partial")
    temp_root = Path(temp_dir) if temp_dir else sample_dir
    # individualID in the prefix keeps concurrent samples apart when they share
    # one --temp-dir, whether from --jobs or from several array tasks on a node.
    temp_prefix = temp_root / f"{sample.individual_id}.sorttmp"
    try:
        if final_bam.is_file() and not overwrite:
            print(f"{sample.individual_id}: {final_bam.name} already present, skipping "
                  f"(pass --overwrite to realign)")
            return _summary_row(sample, final_bam, "skipped")

        fwd, rev = sample.trimmed_wgs_pair(qc_output, paired_tokens)
        sample_dir.mkdir(parents=True, exist_ok=True)
        temp_root.mkdir(parents=True, exist_ok=True)
        _clear_sort_temps(temp_prefix)

        bwa_rc, sort_rc = _bwa_sort(ref, fwd, rev, partial_bam, temp_prefix,
                                    threads, sort_threads, sort_memory)
        if bwa_rc != 0 or sort_rc != 0:
            _clear_sort_temps(temp_prefix)
            if partial_bam.is_file():
                partial_bam.unlink()
            return _summary_row(
                    sample, None, "FAILED",
                    error=f"bwa mem exited {bwa_rc}, samtools sort exited {sort_rc}",
                    )
        if not partial_bam.is_file():
            return _summary_row(sample, None, "FAILED",
                                error="pipeline exited 0 but no BAM was written")

        # Rename only after both halves of the pipe succeeded, so a killed job
        # never leaves a truncated BAM that skip-if-exists would accept.
        os.replace(partial_bam, final_bam)
        index_result = subprocess.run(["samtools", "index", str(final_bam)])
        if index_result.returncode != 0:
            return _summary_row(sample, final_bam, "FAILED",
                                error=f"samtools index exited {index_result.returncode}")
        return _summary_row(sample, final_bam, "aligned")
    except Exception as exc:
        return _summary_row(sample, None, "FAILED", error=f"{type(exc).__name__}: {exc}")


def run_alignment(args, samples, ref, total_samples):
    """
    Across --jobs samples.
    """
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
                pool.submit(_align_one, sample, args.qc_output, output_dir, ref,
                            args.temp_dir, args.paired, args.threads,
                            args.sort_threads, args.sort_memory, args.overwrite)
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
    print(f"\nWGS alignment: {len(aligned)} aligned, {len(skipped)} skipped, "
          f"{len(failed)} failed")
    print(f"Wrote summary: {summary_csv}")
    if args.sample_index is not None:
        _merge_summary_parts(output_dir, total_samples)
    for row in failed:
        print(f"  FAILED {row['individualID']}: {row['error']}", file=sys.stderr)
    if failed:
        raise SystemExit(
                f"wgseqpipe align: {len(failed)} sample(s) failed -- see {summary_csv}"
                )
    return rows


def run(args):
    if args.hpc:
        _submit_hpc(args)
        return
    meta = load_metadata(args.metadata)
    samples = _select_samples(meta, args.sample_index)
    fasta = ensure_unzipped_fasta(args.ref, work_dir=args.work_dir)
    ref = _resolve_bwa_index(fasta, args.build_index)
    print(f"Reference: {ref}")
    print(f"BWA-mem: {len(samples)} of {len(meta.samples)} sample(s) "
          f"({', '.join(s.individual_id for s in samples)}), {args.jobs} at a time, "
          f"{args.threads} bwa threads + {args.sort_threads} sort threads each")
    print(f"samtools sort temp root: {args.temp_dir or '<each sample output directory>'}")
    run_alignment(args, samples, ref, total_samples=len(meta.samples))


def _submit_hpc(args):
    command_parts = [
        # -u: stdout is block-buffered, so a wall-clock kill otherwise discards
        # every print in this module -- including the line saying what --ref
        # resolved to and whether its index was found.
        "python", "-u", "-m", "tarva.cli", "wgseqpipe", "align",
        "--metadata", str(args.metadata),
        "--qc-output", str(args.qc_output),
        "--output", str(args.output),
        "--ref", str(args.ref),
        "--paired", ",".join(args.paired),
        "--threads", str(args.threads),
        "--sort-threads", str(args.sort_threads),
        "--sort-memory", str(args.sort_memory),
        "--jobs", str(args.jobs),
        ]
    if args.work_dir:
        command_parts += ["--work-dir", str(args.work_dir)]
    if args.temp_dir:
        command_parts += ["--temp-dir", str(args.temp_dir)]
    if args.build_index:
        command_parts.append("--build-index")
    if args.overwrite:
        command_parts.append("--overwrite")

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    sbatch_kwargs.setdefault("cpus_per_task", (args.threads + args.sort_threads) * args.jobs)

    log_suffix = ""
    if args.array:
        if args.jobs > 1:
            print(f"WARNING: --array already gives each sample its own task, so --jobs "
                  f"{args.jobs} multiplies every task's core and memory demand. Use "
                  f"--jobs 1 with --array.", file=sys.stderr)
        n_samples = len(load_metadata(args.metadata).samples)
        array_spec = f"0-{n_samples - 1}"
        if args.array_concurrency:
            array_spec += f"%{args.array_concurrency}"
        extra_directives = dict(extra_directives)
        extra_directives.setdefault("array", array_spec)
        # Left unexpanded on purpose: the shell inside the job substitutes it.
        command_parts += ["--sample-index", "$SLURM_ARRAY_TASK_ID"]
        log_suffix = "_%a"
        print(f"Submitting a job array of {n_samples} task(s): --array={array_spec}")

    log_dir = Path(args.output) / "logs"
    script_path = write_sbatch_script(
            job_name="tarva_wgs_align",
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / "tarva_wgs_align.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            log_suffix=log_suffix,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
