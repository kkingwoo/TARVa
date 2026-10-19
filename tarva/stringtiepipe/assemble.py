"""
stringtiepipe assemble: the 3-round StringTie process

    `tarva stringtiepipe assemble --metadata <meta.csv> --aligned <firstpass_dir>
        --output <assembly_dir> --gtf <reference_annotation.gtf[.gz]>
        --strandedness rf|fr|unstranded
        [--stringtie-bin PATH] [--steps STEPS] [--min-transcript-len N]
        [--threads N] [--jobs N] [--overwrite] [--work-dir DIR]
        [--hpc] [--conda-env NAME] [--modules M1,M2] [--directives-dict "{...}"]`

Three rounds:
    1.) round1 - StringTie per sample against the reference annotation
        (`-t -p -G -A -C -B` plus the --strandedness flag), writing
        <id>_round1.gtf
    2.) merge  - every sample's round-1 GTF merged into one updated
        transcriptome (`--merge -G`), writing merged.gtf. This is a hard
        barrier: it consumes all samples at once, so round1 must be complete
        across the whole cohort before it runs.
    3.) round3 - StringTie per sample again, this time against merged.gtf and
        in estimation-only mode (`-m 50 -t -p -G -A -C -B -e` plus the
        --strandedness flag), writing
        <id>_round3.gtf plus the Ballgown .ctab tables downstream reads.

Input:  the first-pass BAMs from `stringtiepipe align`, and the reference
        annotation as a .gff -- StringTie and gffcompare share a parser that
        rejects an NCBI .gtf. --strandedness is REQUIRED and has no default:
        the wrong choice produces wrong assemblies with no error, so fill it
        from the library prep.
Output: <output>/<group>/<individualID>/ per sample, with merged.gtf at the
        <output> root. Each sample gets its own directory so the Ballgown
        .ctab files `-B` emits cannot collide; those tables are what
        `dbpipe load-expression` reads.

StringTie may live in its own conda environment, whose Python cannot run this
CLI. The job runs under --conda-env and invokes StringTie by absolute path via
--stringtie-bin rather than switching environments.
"""

import csv
import subprocess
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tarva.metadata import load_metadata
from tarva.reference import ensure_unzipped
from tarva.stringtiepipe.align import FIRSTPASS_BAM_SUFFIX
from tarva.slurm_utils import parse_directives_dict, split_directives, submit_job, write_sbatch_script

DEFAULT_MIN_TRANSCRIPT_LEN = 50
DEFAULT_STRINGTIE_BIN = "stringtie"
STEP_NAMES = ("round1", "merge", "round3")
MERGED_GTF_NAME = "merged.gtf"
STRANDEDNESS_FLAGS = {"rf": ["--rf"], "fr": ["--fr"], "unstranded": []}

SUMMARY_FIELDS = [
    "individualID", "group", "round", "status", "gtf", "abundance",
    "cov_gtf", "transcripts", "error",
    ]


def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--aligned", required=True,
                        help="The --output directory `stringtiepipe align` wrote; first-pass "
                             "BAMs are read from <aligned>/<group>/<individualID>/")
    parser.add_argument("--output", required=True,
                        help="Directory for assembly output, as <output>/<group>/<individualID>/, "
                             f"with {MERGED_GTF_NAME} at the root")
    parser.add_argument("--gtf", required=True,
                        help="Reference annotation for round 1 and the merge; decompressed "
                             "first if it ends in .gz. Pass the GFF3 here, NOT the GTF given "
                             "to `stringtiepipe align` -- StringTie rejects the NCBI GTF with "
                             "'Error: no valid ID found for GFF record' and exits within "
                             "seconds for every sample, while STAR accepts it fine. NCBI "
                             "names the GFF3 file *_genomic.gff.gz, not .gff3")
    parser.add_argument("--strandedness", required=True, choices=("rf", "fr", "unstranded"),
                        help="Library strandedness for rounds 1 and 3. 'rf' passes StringTie "
                             "--rf (fr-firststrand, e.g. dUTP protocols), 'fr' passes --fr "
                             "(fr-secondstrand), 'unstranded' passes neither. No default: the "
                             "wrong choice produces wrong assemblies without any error, so it "
                             "must be stated explicitly. The merge does not use this option")
    parser.add_argument("--stringtie-bin", default=DEFAULT_STRINGTIE_BIN,
                        help="StringTie executable (default: 'stringtie' from PATH). Pass an "
                             "absolute path when StringTie lives in a different conda env "
                             "than --conda-env, e.g. ~/.conda/envs/stringtie/bin/stringtie")
    parser.add_argument("--steps", default="all",
                        help="Comma-separated subset of rounds to run: "
                             + ", ".join(STEP_NAMES)
                             + " (default: all). Skipping a round requires its output to "
                               "already be present in --output, since each round consumes "
                               "the previous one's")
    parser.add_argument("--min-transcript-len", type=int, default=DEFAULT_MIN_TRANSCRIPT_LEN,
                        help=f"StringTie -m for round 3 (default: {DEFAULT_MIN_TRANSCRIPT_LEN}); "
                             "round 1 does not set -m, matching the methods")
    parser.add_argument("--work-dir",
                        help="Where to write the decompressed GTF if --gtf is gzipped "
                             "(default: alongside --gtf)")
    parser.add_argument("--threads", type=int, default=8,
                        help="StringTie -p per sample (default: 8)")
    parser.add_argument("--jobs", type=int, default=2,
                        help="Samples to assemble concurrently (default: 2). StringTie is far "
                             "lighter on memory than STAR -- there is no genome index to load "
                             "-- so concurrency is cheap here")
    parser.add_argument("--overwrite", action="store_true",
                        help="Re-run rounds for samples that already have output")
    parser.add_argument("--hpc", action="store_true",
                        help="Submit as a SLURM job instead of running locally")
    parser.add_argument("--conda-env", default="bio_qc",
                        help="Conda environment activated in the generated SLURM script "
                             "(default: bio_qc -- runs the CLI itself)")
    parser.add_argument("--modules", default="",
                        help="Comma-separated HPC modules to load (default: none -- StringTie "
                             "comes from --stringtie-bin, not a module)")
    parser.add_argument("--directives-dict", dest="directives_dict",
                        type=parse_directives_dict, default={}, metavar="DICT",
                        help="SLURM sbatch directives as a Python dict literal, only used with "
                             "--hpc. Example: \"{'partition': 'Orion', 'mem': '32gb', "
                             "'time': '12:00:00', 'cpus_per_task': 8}\"")


def _parse_steps(raw):
    """Same contract as qc.py's --steps: reject typos rather than silently no-op."""
    if str(raw).strip().lower() == "all":
        return set(STEP_NAMES)
    steps = {s.strip().lower() for s in str(raw).split(",") if s.strip()}
    unknown = sorted(steps - set(STEP_NAMES))
    if unknown:
        raise SystemExit(
                f"stringtiepipe assemble: unknown --steps value(s): {unknown}. "
                f"Valid rounds are {list(STEP_NAMES)}, or 'all'."
                )
    if not steps:
        raise SystemExit("stringtiepipe assemble: --steps was empty -- give at least one round.")
    return steps


def _sample_dir(root, sample):
    return Path(root) / sample.group / sample.individual_id


def _firstpass_bam(sample, aligned_root):
    return (_sample_dir(aligned_root, sample)
            / f"{sample.individual_id}{FIRSTPASS_BAM_SUFFIX}")


def _round_paths(sample, output_dir, round_label):
    sample_dir = _sample_dir(output_dir, sample)
    return {
        "dir": sample_dir,
        "gtf": sample_dir / f"{sample.individual_id}_{round_label}.gtf",
        "cov_gtf": sample_dir / f"{sample.individual_id}_{round_label}_cov.gtf",
        "abundance": sample_dir / f"{sample.individual_id}_{round_label}.tab",
        }


def _count_transcripts(gtf):
    """Transcript records in a StringTie GTF -- a cheap sanity number per sample."""
    gtf = Path(gtf)
    if not gtf.is_file():
        return ""
    count = 0
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) > 2 and fields[2] == "transcript":
                count += 1
    return count


def _summary_row(sample, round_label, status, paths, error=""):
    return {
        "individualID": sample.individual_id,
        "group": sample.group,
        "round": round_label,
        "status": status,
        "gtf": str(paths["gtf"]) if status != "FAILED" else "",
        "abundance": str(paths["abundance"]) if status != "FAILED" else "",
        "cov_gtf": str(paths["cov_gtf"]) if status != "FAILED" else "",
        "transcripts": _count_transcripts(paths["gtf"]) if status != "FAILED" else "",
        "error": error,
        }


def _assemble_one(sample, aligned_root, output_dir, stringtie_bin, annotation, strandedness,
                  threads, round_label, min_transcript_len=None, estimate_only=False,
                  overwrite=False):
    """
    One sample, one round. Failures become FAILED rows rather than exceptions,
    so one bad sample doesn't abort a pool with hours of work in flight.
    """
    paths = _round_paths(sample, output_dir, round_label)
    try:
        if paths["gtf"].is_file() and not overwrite:
            print(f"{sample.individual_id} [{round_label}]: {paths['gtf'].name} already "
                  f"present, skipping (pass --overwrite to redo)")
            return _summary_row(sample, round_label, "skipped", paths)

        bam = _firstpass_bam(sample, aligned_root)
        if not bam.is_file():
            return _summary_row(
                    sample, round_label, "FAILED", paths,
                    error=f"first-pass BAM not found: {bam}",
                    )
        paths["dir"].mkdir(parents=True, exist_ok=True)

        cmd = [
            str(stringtie_bin), str(bam),
            "-o", str(paths["gtf"]),
            *STRANDEDNESS_FLAGS[strandedness],
            "-t",
            "-p", str(threads),
            "-G", str(annotation),
            "-A", str(paths["abundance"]),
            "-C", str(paths["cov_gtf"]),
            "-B",
            ]
        if min_transcript_len is not None:
            cmd += ["-m", str(min_transcript_len)]
        if estimate_only:
            cmd.append("-e")

        print(f"{sample.individual_id} [{round_label}]: " + " ".join(cmd))
        result = subprocess.run(cmd)
        if result.returncode != 0:
            return _summary_row(
                    sample, round_label, "FAILED", paths,
                    error=f"stringtie exited {result.returncode}",
                    )
        if not paths["gtf"].is_file():
            return _summary_row(
                    sample, round_label, "FAILED", paths,
                    error=f"stringtie exited 0 but {paths['gtf'].name} is missing",
                    )
        return _summary_row(sample, round_label, "assembled", paths)
    except Exception as exc:
        return _summary_row(sample, round_label, "FAILED", paths,
                            error=f"{type(exc).__name__}: {exc}")


def _run_round(args, meta, round_label, annotation, min_transcript_len, estimate_only):
    """
    Fan a round out across --jobs samples, flushing the summary per sample so a
    wall-clock kill still leaves a record of what finished.
    """
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = output_dir / f"{round_label}_summary.csv"
    rows = []

    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        fh.flush()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [
                pool.submit(_assemble_one, sample, args.aligned, output_dir,
                            args.stringtie_bin, annotation, args.strandedness,
                            args.threads, round_label,
                            min_transcript_len, estimate_only, args.overwrite)
                for sample in meta.samples
                ]
            for future in as_completed(futures):
                row = future.result()
                writer.writerow(row)
                fh.flush()
                rows.append(row)

    failed = [r for r in rows if r["status"] == "FAILED"]
    done = [r for r in rows if r["status"] in ("assembled", "skipped")]
    print(f"{round_label}: {len(done)} complete, {len(failed)} failed -- {summary_csv}")
    for row in failed:
        print(f"  FAILED {row['individualID']}: {row['error']}", file=sys.stderr)
    if failed:
        raise SystemExit(
                f"stringtiepipe assemble: {len(failed)} sample(s) failed in {round_label} "
                f"-- see {summary_csv}. The merge needs every sample, so fix these before "
                f"continuing."
                )
    return rows


def _collect_round1_gtfs(meta, output_dir):
    """
    The round-1 GTFs the merge consumes, built from expected per-sample paths.

    Deliberately not a directory glob: the `-C` coverage output is also a
    .gtf, so a glob would merge it in. Constructing the list from the sample
    roster makes that impossible and also catches a missing sample rather than
    quietly merging fewer.
    """
    gtfs, missing = [], []
    for sample in meta.samples:
        gtf = _round_paths(sample, output_dir, "round1")["gtf"]
        (gtfs if gtf.is_file() else missing).append(gtf)
    if missing:
        raise SystemExit(
                "stringtiepipe assemble: round-1 GTF missing for "
                f"{len(missing)} sample(s): {[str(m) for m in missing]}. Run --steps round1 "
                "first -- merging a partial cohort would silently change the transcriptome."
                )
    return gtfs


def run_merge(stringtie_bin, annotation, round1_gtfs, output_dir):
    """Round 2: one call consuming every sample's round-1 GTF."""
    merged = Path(output_dir) / MERGED_GTF_NAME
    cmd = ([str(stringtie_bin), "--merge", "-G", str(annotation), "-o", str(merged)]
           + [str(g) for g in round1_gtfs])
    print(f"merge ({len(round1_gtfs)} sample GTFs): " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    if not merged.is_file():
        raise SystemExit(f"stringtiepipe assemble: merge exited 0 but {merged} is missing")
    print(f"Wrote merged transcriptome: {merged} "
          f"({_count_transcripts(merged)} transcripts)")
    return merged


def run(args):
    if args.hpc:
        _submit_hpc(args)
        return

    steps = _parse_steps(args.steps)
    meta = load_metadata(args.metadata)
    output_dir = Path(args.output)
    merged_gtf = output_dir / MERGED_GTF_NAME
    reference_gtf = ensure_unzipped(args.gtf, work_dir=args.work_dir,
                                    description="annotation GTF")
    print("Running StringTie round(s): " + ", ".join(s for s in STEP_NAMES if s in steps))

    if "round1" in steps:
        _run_round(args, meta, "round1", reference_gtf,
                   min_transcript_len=None, estimate_only=False)

    if "merge" in steps:
        run_merge(args.stringtie_bin, reference_gtf,
                  _collect_round1_gtfs(meta, output_dir), output_dir)
    elif "round3" in steps and not merged_gtf.is_file():
        raise SystemExit(
                f"stringtiepipe assemble: --steps includes 'round3' but {merged_gtf} does not "
                f"exist. Round 3 estimates against the merged transcriptome, so run 'merge' "
                f"at least once first."
                )

    if "round3" in steps:
        _run_round(args, meta, "round3", merged_gtf,
                   min_transcript_len=args.min_transcript_len, estimate_only=True)

    print(f"Assembly complete: {output_dir}")


def _submit_hpc(args):
    command_parts = [
        # -u: unbuffered, so a job killed at the wall clock still leaves its
        # progress in the .o rather than discarding the buffer.
        "python", "-u", "-m", "tarva.cli", "stringtiepipe", "assemble",
        "--metadata", str(args.metadata),
        "--aligned", str(args.aligned),
        "--output", str(args.output),
        "--gtf", str(args.gtf),
        "--strandedness", str(args.strandedness),
        "--stringtie-bin", str(args.stringtie_bin),
        "--steps", str(args.steps),
        "--min-transcript-len", str(args.min_transcript_len),
        "--threads", str(args.threads),
        "--jobs", str(args.jobs),
        ]
    if args.work_dir:
        command_parts += ["--work-dir", str(args.work_dir)]
    if args.overwrite:
        command_parts.append("--overwrite")

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    sbatch_kwargs.setdefault("cpus_per_task", args.threads * args.jobs)
    log_dir = Path(args.output) / "logs"
    script_path = write_sbatch_script(
            job_name="tarva_stringtie_assemble",
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / "tarva_stringtie_assemble.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
