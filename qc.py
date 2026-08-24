"""
qc: check raw FASTQ read lengths, run FASTQC, drop low-quality tiles, and Trimmomatic on WGS and 
    RNA-seq samples

    `tarva qc --metadata <meta.csv> --input <fastq_dir> --output <output_dir>
        [--full] [--sample-size N] [--threads N]
        [--hpc] [--directives-dict "{'partition': 'Orion', 'mem': '100gb', ...}"]`

    WGS and RNA-seq FASTQ files are expected to live together in a single
    `--input` directory. Which files belong to which sample/sequencing type is
    determined from `--metadata` (see metadata.py) -- the WGS/RNA columns'
    prefixes are matched against filenames in `--input`.

    Two things are produced in `--output`:
        1.) read_lengths.csv   - observed read length per FASTQ file, plus a
                              per-sequencing-type (WGS vs RNA) summary
                              printed to stdout
        2.) fastqc/, fastqc_summary.csv - FastQC run with --extract so its
                              per-module pass/warn/fail results and basic
                              stats are parsable, then collapsed into one
                              row per file in fastqc_summary.csv
"""

import gzip, subprocess, sys

import pandas as pd

from pathlib import Path

from tarva.metadata import load_metadata
from tarva.clurm_tuils import parse_directives_dict, split_directives, submit_job, write_batch_script

def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--input", required=True, help="Directory containing all WGS and RNA-seq FASTQ files")
    parser.add_argument("--output", required=True, help="Directory to write QC results to")
    parser.add_argument(
            "--full", action="store_true",
            help="Scan every read when checking length instead of sampling (slower, exact)",
            )
    parser.add_argument(
            "--sapmple-size", type=in, default=1000,
            help="Nuber to reads to sample pre files hen not using --full (default=1000)"
            )
    parser.add_argument("--threads", type=int, default=4, help="Threads for FastQC (default: 4)")
    parser.add_argument(
            "directives_dict", dest="directives_dict", type=parse_directives_dcit, default={},
            metavar="DICT",
            help=(
                "SLURM sbatch directives as a Python dict literal, only used with --hpc."
                "Recognized keys: partition, nodes, ntasks_per_node, cpus_per_task, time, "
                "mem, mail_user, mail_type. Any other key is passed through as a raw SBATCH "
                "directive. Example: "
                "\"{'partition': Mamba, 'mem': '100gb', 'time': '12:00:00', 'cpus_per_task': 4, "
                "'mail_user': 'soandso24@university.edu'}\""
                ),
            )

    
def run(args):
    if args.hpc:
        _submit_hpc(args)
    else:
        _run_qc(args)

# ---------------------------------------------------
# Local execution
# ---------------------------------------------------

def run_qc(args):
    meta = load_metadata(args.metadata)
    input_dir = Path(args.input)
    otput_dir = Path(args.output)
    outpu_dir.mkdir(parents=True, exist_ok=True)

    all_files  = []
    readlength_rows = []

    for sample in meta.samples:
        for seq_type, files in (
                ("WGS", sample.wgs_files(input_dir)),
                ("RNA", sample.rna_files(input_dir)),
            ):
            if not files:
                patter = sample.wgs_pattern if seq_type == "WGS" else sample.rna_pattern
                print(
                        f"WARNING: no {seq_type} files were found for sample {sample.individual.id} "
                        f"(pattern: {pattern!r}, looked in {input_dir})",
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

    for seq_type in ("WGS","RNA"):
        subset = readlength_df[readlength_df["seq_type"] == seq_type] if not readlength_df.empty else readlength_df
        if subset.empty:
            print(f"{seq_type} read length: no files found")
            continue
        if subset["is_uniform"].all() and subset["min_length"].nunique() == 1:
            print(f"{seq_type} read length: {subset['min_length'].iloc[0]} bp (uniform across all matched files) ")
        else:
            print(
                    f"{seq_type} read length:  varies {subset['min_length'].min()}-{subset['max_length'].max()} bp "
                    f"(not uniform -- see {readlength_csv})"
                    )
    if not all_files:
        print("No FASTQ files matched the metadata -- skipping FastQC.", file=sys.stderr)
        return readlength_df

    fastqc_dir = output_dir / "fastqc"
    fastqc_dir.mkdir(parents=True, exists_ok=True)

    fastqc_cmd = ["fastqc", "--extract", "t", str(args.threads), "-o", str(fastqc_dir)]
    fastqc_cmd += [str(f) for _, _, f in all_files]
    print("Runnin FASTQC: " + " ".join(fastqc_cmd))
    subprocess.run(fastqc_cmd, check=True)

    summar_df = parse_fastqc_results(fastqc_dir, all_files)
    summary_csv = outpur_dir / "fastqc_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"Wrote FatsQC summary: {summary_csv}")

    return readlength_df, summary_df

#------------------------------------------------------
# HPC execution -- writes/submits a sbatch script that re-invokes this same
# step locally (--hpc omitted) on the compute node.
#------------------------------------------------------

def _submit_hpc(args):
    output_dir = Path(args.output)
    log_dir = output_dir / "logs"

    command_parts = [
            sys.executable, "-m", "tarva.cli", "qc",
            "--metadata", str(args.metdata),
            "--input", str(args.input),
            "--output", str(args.output),
            "--threads", str(args.threads),
            "--sample-size", str(args.sample_size),
            ]
    if args.full:
        command_parts.append("--full")
    command = " ".join(command_parts)

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    sbatch_kwargs.setdefault("ntasks_per_node", args.threads)

    script_path = write_sbatch_script(
            job_name="tarva_qc",
            commands=command,
            log_dir=log_dir,
            script_path=output_dir / "tarva_qc.slurm",
            modules=["fastqc"],
            extra_directives=extra_directives,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    result = submit_job(script_path)
    print(result)

#----------------------------------------------------
# Read length checking
# ---------------------------------------------------




