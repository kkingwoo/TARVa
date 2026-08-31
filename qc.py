"""
qc: check raw FASTQ read lengths, run FastQC, and trim reads with Trimmomatic
    for WGS and RNA-seq samples

    `tarva qc --metadata <meta.csv> --input <fastq_dir> --output <output_dir>
        [--full] [--sample-size N] [--threads N] [--jobs N]
        [--minlen-fraction F] [--list-adapters]
        [--hpc] [--conda-env NAME] [--directives-dict "{'partition': ..., ...}"]`

    FASTQ files are expected to live in a directory tree under `--input`,
    produced by `tarva sort-input`: <input>/WGS/<group>/<adapter_stem>/... and
    <input>/RNA/<group>/<adapter_stem>/... (see metadata.py for details).
    Which files belong to which sample within that tree is determined from
    `--metadata` -- the WGS/RNA columns' prefixes are matched against
    filenames in the resolved subfolder.

    Steps:
        1.) Observed read length is checked per FASTQ file (read_lengths.csv).
        2.) FastQC is run per (seq_type, group, adapter) group, writing into
            the mirrored tree <output>/fastqc/WGS|RNA/<group>/<adapter_stem>/
            (fastqc_summary.csv summarizes across all of them). Up to --jobs
            groups are run concurrently, each FastQC invocation itself using
            --threads internally.
        3.) Trimmomatic is run per paired sample, using the observed read
            length to set MINLEN (round(--minlen-fraction * read_length)).
            The adapter FASTA used for each sample/seq-type comes from its
            wgs_barcode/rna_barcode metadata value. Output mirrors the same
            tree: <output>/trimmed/WGS|RNA/<group>/<adapter_stem>/. Up to
            --jobs samples are trimmed concurrently, each Trimmomatic
            invocation itself using --threads internally.

    Re-run as needed until reads look clean -- each step's output directory is
    overwritten on each run.
"""

import gzip, subprocess, sys

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from tarva.metadata import load_metadata
from tarva.adapters import print_available_adapters, resolve_adapter
from tarva.slurm_utils import parse_directives_dict, split_directives, submit_job, write_sbatch_script

DEFAULT_MINLEN_FRACTION = 0.7


def add_arguments(parser):
    parser.add_argument(
            "--metadata",
            help="Path to metadata CSV (see metadata.py) -- required unless --list-adapters is used",
            )
    parser.add_argument(
            "--input",
            help=(
                "Root directory of the sorted FASTQ tree produced by `tarva sort-input` "
                "(<input>/WGS/<group>/<adapter>/... and <input>/RNA/<group>/<adapter>/...) "
                "-- required unless --list-adapters is used"
                ),
            )
    parser.add_argument(
            "--output",
            help="Directory to write QC results to -- required unless --list-adapters is used",
            )
    parser.add_argument(
            "--full", action="store_true",
            help="Scan every read when checking length instead of sampling (slower, exact)",
            )
    parser.add_argument(
            "--sample-size", type=int, default=1000,
            help="Number of reads to sample per file when not using --full (default: 1000)",
            )
    parser.add_argument("--threads", type=int, default=4, help="Threads for FastQC/Trimmomatic (default: 4)")
    parser.add_argument(
            "--jobs", type=int, default=1,
            help=(
                "Number of FastQC groups / samples' Trimmomatic runs to execute "
                "concurrently (default: 1). Each job uses --threads internally, so total "
                "CPU need is roughly --jobs x --threads -- size --cpus-per-task in "
                "--directives-dict accordingly when using --hpc."
                ),
            )
    parser.add_argument(
            "--minlen-fraction", type=float, default=DEFAULT_MINLEN_FRACTION,
            help=(
                "Trimmomatic MINLEN is set to round(fraction * observed read length) "
                f"per file (default: {DEFAULT_MINLEN_FRACTION})"
                ),
            )
    parser.add_argument(
            "--list-adapters", action="store_true",
            help="List the bundled adapter FASTA files, then exit",
            )
    parser.add_argument(
            "--hpc", action="store_true",
            help="Submit this step as a SLURM job instead of running locally",
            )
    parser.add_argument(
            "--conda-env", default="bio_qc",
            help="Conda environment to activate in the generated SLURM script (default: bio_qc)",
            )
    parser.add_argument(
            "--directives-dict", dest="directives_dict", type=parse_directives_dict, default={},
            metavar="DICT",
            help=(
                "SLURM sbatch directives as a Python dict literal, only used with --hpc. "
                "Recognized keys: partition, nodes, ntasks_per_node, cpus_per_task, time, "
                "mem, mail_user, mail_type. Any other key is passed through as a raw SBATCH "
                "directive. Example: "
                "\"{'partition': 'Mamba', 'mem': '100gb', 'time': '12:00:00', 'cpus_per_task': 4, "
                "'mail_user': 'soandso24@university.edu'}\""
                ),
            )


def run(args):
    if args.list_adapters:
        print_available_adapters()
        return
    missing = [name for name in ("metadata", "input", "output") if getattr(args, name) is None]
    if missing:
        raise SystemExit(
                "tarva qc: the following arguments are required: "
                + ", ".join(f"--{name}" for name in missing)
                )
    if args.hpc:
        _submit_hpc(args)
    else:
        run_qc(args)

# ---------------------------------------------------
# Local execution
# ---------------------------------------------------

def _write_output_location_note(output_dir):
    note_path = Path.cwd() / "tarva_qc_output_location.txt"
    lines = [
            f"TARVa QC results and logs: {output_dir}",
            "  read_lengths.csv          - per-file read length summary",
            "  fastqc/                   - FastQC output, organized WGS|RNA/<group>/<adapter>/",
            "  fastqc_summary.csv        - collapsed FastQC summary, one row per file",
            "  trimmed/                  - Trimmomatic output, organized WGS|RNA/<group>/<adapter>/",
            "  trimmomatic_summary.csv   - one row per sample/seq-type trimmed",
            "  logs/                     - SLURM job logs (only present when run with --hpc)",
            ]
    note_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote output location summary: {note_path}")

def run_qc(args):
    meta = load_metadata(args.metadata)
    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_output_location_note(output_dir)

    all_files = []
    readlength_rows = []

    for sample in meta.samples:
        for seq_type, files in (
                ("WGS", sample.wgs_files(input_dir)),
                ("RNA", sample.rna_files(input_dir)),
            ):
            if not files:
                pattern = sample.wgs_pattern if seq_type == "WGS" else sample.rna_pattern
                seq_dir = sample.wgs_dir(input_dir) if seq_type == "WGS" else sample.rna_dir(input_dir)
                print(
                        f"WARNING: no {seq_type} files were found for sample {sample.individual_id} "
                        f"(pattern: {pattern!r}, looked in {seq_dir})",
                        file=sys.stderr,
                        )
                continue
            for f in files:
                all_files.append((sample, seq_type, f))
                stats = check_read_length(f, full=args.full, sample_size=args.sample_size)
                readlength_rows.append({
                    "individualID": sample.individual_id,
                    "group": sample.group,
                    "seq_type": seq_type,
                    "file": f.name,
                    "n_reads_checked": stats["n_reads"],
                    "min_length": stats["min"],
                    "max_length": stats["max"],
                    "is_uniform": stats["min"] == stats["max"],
                    })

    readlength_df = pd.DataFrame(readlength_rows)
    readlength_csv = output_dir / "read_lengths.csv"
    readlength_df.to_csv(readlength_csv, index=False)
    print(f"Wrote read length summary: {readlength_csv}")

    for seq_type in ("WGS", "RNA"):
        subset = readlength_df[readlength_df["seq_type"] == seq_type] if not readlength_df.empty else readlength_df
        if subset.empty:
            print(f"{seq_type} read length: no files found")
            continue
        if subset["is_uniform"].all() and subset["min_length"].nunique() == 1:
            print(f"{seq_type} read length: {subset['min_length'].iloc[0]} bp (uniform across all matched files)")
        else:
            print(
                    f"{seq_type} read length: varies {subset['min_length'].min()}-{subset['max_length'].max()} bp "
                    f"(not uniform -- see {readlength_csv})"
                    )

    if not all_files:
        print("No FASTQ files matched the metadata -- skipping FastQC and Trimmomatic.", file=sys.stderr)
        return readlength_df

    fastqc_groups = defaultdict(list)
    for sample, seq_type, f in all_files:
        barcode = sample.wgs_barcode if seq_type == "WGS" else sample.rna_barcode
        adapter_stem = resolve_adapter(barcode).stem
        fastqc_groups[(seq_type, sample.group, adapter_stem)].append((sample, seq_type, f))

    def _run_fastqc_group(group_key):
        seq_type, group, adapter_stem = group_key
        group_dir = _fastqc_group_dir(output_dir, seq_type, group, adapter_stem)
        group_dir.mkdir(parents=True, exist_ok=True)
        cmd = ["fastqc", "--extract", "-t", str(args.threads), "-o", str(group_dir)]
        cmd += [str(f) for _, _, f in fastqc_groups[group_key]]
        print("Running FastQC: " + " ".join(cmd))
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(
                    f"WARNING: FastQC reported one or more failures for group {group_key} "
                    f"(exit code {result.returncode}) -- see the FastQC output above for which "
                    f"file(s) failed. Continuing with whatever files it did process successfully.",
                    file=sys.stderr,
                    )
        return group_key

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(_run_fastqc_group, key) for key in fastqc_groups]
        for future in as_completed(futures):
            future.result()

    summary_df = parse_fastqc_results(output_dir, all_files)
    summary_csv = output_dir / "fastqc_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"Wrote FastQC summary: {summary_csv}")

    trimmed_df = run_trimmomatic(
            meta, input_dir, output_dir, readlength_df,
            minlen_fraction=args.minlen_fraction,
            threads=args.threads,
            jobs=args.jobs,
            )
    trimmed_csv = output_dir / "trimmomatic_summary.csv"
    trimmed_df.to_csv(trimmed_csv, index=False)
    print(f"Wrote Trimmomatic summary: {trimmed_csv}")

    return readlength_df, summary_df, trimmed_df

# ---------------------------------------------------
# Read length checking
# ---------------------------------------------------

def _open_fastq(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")

def check_read_length(path, full=False, sample_size=1000):
    """
    Returns {"n_reads": int, "min": int, "max": int} for the sequence lines
    (every 2nd of each 4-line FASTQ record) in `path`. Reads are counted
    from the start of the file; when `full` is False, only the first
    `sample_size` reads are inspected.
    """
    min_len = None
    max_len = None
    n_reads = 0
    with _open_fastq(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            length = len(line.rstrip("\n"))
            min_len = length if min_len is None else min(min_len, length)
            max_len = length if max_len is None else max(max_len, length)
            n_reads += 1
            if not full and n_reads >= sample_size:
                break
    return {"n_reads": n_reads, "min": min_len, "max": max_len}

# ---------------------------------------------------
# FastQC summary parsing
# ---------------------------------------------------

def _fastqc_group_dir(output_dir, seq_type, group, adapter_stem):
    return output_dir / "fastqc" / seq_type / group / adapter_stem

def _fastqc_extract_dir(output_dir, sample, seq_type, fastq_path):
    barcode = sample.wgs_barcode if seq_type == "WGS" else sample.rna_barcode
    adapter_stem = resolve_adapter(barcode).stem
    group_dir = _fastqc_group_dir(output_dir, seq_type, sample.group, adapter_stem)

    name = Path(fastq_path).name
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    extract_dir = group_dir / f"{name}_fastqc"
    if not extract_dir.is_dir():
        raise FileNotFoundError(
                f"Expected FastQC output directory not found: {extract_dir} "
                f"(looked for extracted results of {fastq_path})"
                )
    return extract_dir

def _parse_basic_statistics(fastqc_data_path):
    stats = {}
    in_block = False
    with open(fastqc_data_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">>Basic Statistics"):
                in_block = True
                continue
            if in_block and line.startswith(">>END_MODULE"):
                break
            if in_block and "\t" in line:
                key, value = line.split("\t", 1)
                stats[key] = value
    return stats

def _parse_module_summary(summary_path):
    counts = {"PASS": 0, "WARN": 0, "FAIL": 0}
    modules = {}
    with open(summary_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            status, module = parts[0], parts[1]
            counts[status] = counts.get(status, 0) + 1
            modules[module] = status
    return counts, modules

def parse_fastqc_results(output_dir, all_files):
    rows = []
    for sample, seq_type, f in all_files:
        try:
            extract_dir = _fastqc_extract_dir(output_dir, sample, seq_type, f)
            basic_stats = _parse_basic_statistics(extract_dir / "fastqc_data.txt")
            counts, modules = _parse_module_summary(extract_dir / "summary.txt")
        except FileNotFoundError:
            print(
                    f"WARNING: no FastQC results found for {f} -- FastQC likely failed to "
                    f"process this file (see the FastQC output above); marking it as failed "
                    f"in the summary instead of skipping it silently.",
                    file=sys.stderr,
                    )
            rows.append({
                "individualID": sample.individual_id,
                "group": sample.group,
                "seq_type": seq_type,
                "file": f.name,
                "total_sequences": None,
                "sequence_length": None,
                "percent_gc": None,
                "n_pass": None,
                "n_warn": None,
                "n_fail": None,
                "failed_modules": "FASTQC_FAILED",
                })
            continue
        rows.append({
            "individualID": sample.individual_id,
            "group": sample.group,
            "seq_type": seq_type,
            "file": f.name,
            "total_sequences": basic_stats.get("Total Sequences"),
            "sequence_length": basic_stats.get("Sequence length"),
            "percent_gc": basic_stats.get("%GC"),
            "n_pass": counts.get("PASS", 0),
            "n_warn": counts.get("WARN", 0),
            "n_fail": counts.get("FAIL", 0),
            "failed_modules": ";".join(m for m, s in modules.items() if s == "FAIL"),
            })
    return pd.DataFrame(rows)

# ---------------------------------------------------
# Trimmomatic
# ---------------------------------------------------

def _minlen_for_file(readlength_df, filename, fraction):
    row = readlength_df[readlength_df["file"] == filename]
    if row.empty:
        raise ValueError(f"No read length recorded for {filename} -- was read length checking skipped?")
    observed = row["max_length"].iloc[0]
    return max(1, round(fraction * observed))

def _trim_one(job):
    sample, seq_type, fwd, rev, barcode, readlength_df, minlen_fraction, threads, trimmed_dir = job
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    minlen = _minlen_for_file(readlength_df, fwd.name, minlen_fraction)
    adapter_path = resolve_adapter(barcode)

    out_fwd_paired = trimmed_dir / f"{fwd.stem}.paired.fastq.gz"
    out_fwd_unpaired = trimmed_dir / f"{fwd.stem}.unpaired.fastq.gz"
    out_rev_paired = trimmed_dir / f"{rev.stem}.paired.fastq.gz"
    out_rev_unpaired = trimmed_dir / f"{rev.stem}.unpaired.fastq.gz"

    cmd = [
        "trimmomatic", "PE", "-threads", str(threads),
        str(fwd), str(rev),
        str(out_fwd_paired), str(out_fwd_unpaired),
        str(out_rev_paired), str(out_rev_unpaired),
        f"ILLUMINACLIP:{adapter_path}:2:30:10",
        "SLIDINGWINDOW:4:20", "LEADING:5", "TRAILING:5", f"MINLEN:{minlen}",
        ]
    print("Running Trimmomatic: " + " ".join(cmd))
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)

    return {
        "individualID": sample.individual_id,
        "group": sample.group,
        "seq_type": seq_type,
        "forward": fwd.name,
        "reverse": rev.name,
        "adapter_used": adapter_path.name,
        "minlen_used": minlen,
        "paired_output": [str(out_fwd_paired), str(out_rev_paired)],
        "log": result.stderr.strip(),
        }

def run_trimmomatic(meta, input_dir, output_dir, readlength_df, minlen_fraction, threads, jobs):
    trim_jobs = []
    for sample in meta.samples:
        for seq_type, files, barcode in (
                ("WGS", sample.wgs_files(input_dir), sample.wgs_barcode),
                ("RNA", sample.rna_files(input_dir), sample.rna_barcode),
            ):
            if len(files) != 2:
                if files:
                    print(
                            f"WARNING: expected exactly 2 paired FASTQ files for "
                            f"{sample.individual_id} ({seq_type}), found {len(files)} -- skipping Trimmomatic",
                            file=sys.stderr,
                            )
                continue
            fwd, rev = files
            adapter_stem = resolve_adapter(barcode).stem
            trimmed_dir = output_dir / "trimmed" / seq_type / sample.group / adapter_stem
            trim_jobs.append(
                    (sample, seq_type, fwd, rev, barcode, readlength_df, minlen_fraction, threads, trimmed_dir)
                    )

    rows = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(_trim_one, job) for job in trim_jobs]
        for future in as_completed(futures):
            rows.append(future.result())

    return pd.DataFrame(rows)

#----------------------------------------------------
# HPC execution -- writes/submits a sbatch script that re-invokes this same
# step locally (--hpc omitted) on the compute node.
#----------------------------------------------------

def _submit_hpc(args):
    output_dir = Path(args.output)
    log_dir = output_dir / "logs"
    _write_output_location_note(output_dir)

    command_parts = [
            "python", "-m", "tarva.cli", "qc",
            "--metadata", str(args.metadata),
            "--input", str(args.input),
            "--output", str(args.output),
            "--threads", str(args.threads),
            "--jobs", str(args.jobs),
            "--sample-size", str(args.sample_size),
            "--minlen-fraction", str(args.minlen_fraction),
            ]
    if args.full:
        command_parts.append("--full")
    command = " ".join(command_parts)

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    sbatch_kwargs.setdefault("ntasks_per_node", args.threads * args.jobs)

    script_path = write_sbatch_script(
            job_name="tarva_qc",
            commands=command,
            log_dir=log_dir,
            script_path=output_dir / "tarva_qc.slurm",
            modules=["anaconda3"],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    result = submit_job(script_path)
    print(result)
