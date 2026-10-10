"""
wgseqpipe pon: build a panel of normals from the control group

    `tarva wgseqpipe pon --metadata <meta.csv> --preprocessed <preproc_dir>
        --output <pon_dir> --ref <genome.fa[.gz]> --control-group <value>
        [--existing-pon <pon.vcf.gz>]
        [--germline-resource <af-only.vcf.gz>] [--intervals <list>]
        [--min-sample-count N] [--steps call,import,create]
        [--max-mnp-distance N] [--java-heap SIZE] [--temp-dir DIR]
        [--jobs N] [--sample-index N] [--overwrite]
        [--hpc] [--array] [--array-concurrency N] [--conda-env NAME]
        [--directives-dict "{...}"]`

Three stages, selectable with --steps, because only the first fans out:

  call    Mutect2 on each control sample
          -> <output>/calls/<individualID>_pon.vcf.gz
  import  GenomicsDBImport of those VCFs  -> <output>/pon_db/
  create  CreateSomaticPanelOfNormals     -> <output>/pon.vcf.gz

With --hpc --array the `call` stage becomes a job array and `import,create` are
submitted as a second job with 'dependency': 'afterok:<array job id>', so the
whole thing is one command.

**Which group is the control is a flag, not a convention.** `--control-group`
is validated against the two values actually present in the metadata, so there
is no hardcoded "Control" string to silently miss. See wgseqpipe/variants.py for
the other half of this: the PON is applied to BOTH groups.

**Building from the control group is the default, not the only option.**
`--existing-pon <vcf>` adopts a panel built elsewhere -- downloaded from a different
source or one from a previous run -- and skips all three stages.
The adopted panel is linked into
<output>/pon.vcf.gz alongside its index, so the command that follows is the same
whichever way the panel was obtained.

"""

import csv
import re
import shutil
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tarva.gatk_common import (DEFAULT_JAVA_HEAP, GatkStepError, run_gatk)
from tarva.metadata import load_metadata
from tarva.reference import ensure_unzipped_fasta
from tarva.slurm_utils import (parse_directives_dict, split_directives, submit_job,
                               write_sbatch_script)

PON_VCF_SUFFIX = "_pon.vcf.gz"
PON_NAME = "pon.vcf.gz"
PON_DB_NAME = "pon_db"
DEFAULT_MIN_SAMPLE_COUNT = 2
ALL_STEPS = ("call", "import", "create")

SUMMARY_FIELDS = ["individualID", "group", "status", "vcf", "error"]


def parse_steps(value):
    steps = [s.strip() for s in str(value).split(",") if s.strip()]
    unknown = [s for s in steps if s not in ALL_STEPS]
    if unknown:
        raise ValueError(f"Unknown --steps value(s) {unknown}; valid: {list(ALL_STEPS)}")
    return steps


def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--preprocessed", required=True,
                        help="The --output directory `wgseqpipe preprocess` wrote; BAMs are "
                             "read from <preprocessed>/<group>/<individualID>/")
    parser.add_argument("--output", required=True, help="Directory for the PON and its inputs")
    parser.add_argument("--ref", required=True,
                        help="Reference genome FASTA; must be the one the BAMs were aligned "
                             "and preprocessed against")
    parser.add_argument("--control-group", required=True, metavar="VALUE",
                        help="Which value of the metadata `group` column supplies the normals. "
                             "Validated against the two values actually present, so a typo is "
                             "an error rather than an empty panel")
    parser.add_argument("--existing-pon", metavar="VCF",
                        help="Adopt an already-built panel of normals -- downloaded from a different "
                             "source or one from an earlier run -- instead of "
                             "building one from --control-group. All three stages are skipped "
                             "and the file is linked to <output>/pon.vcf.gz with its index, so "
                             "`wgseqpipe variants --panel-of-normals` takes the same path "
                             "either way. Must be indexed (.tbi or .idx)")
    parser.add_argument("--bam-suffix", default="_bqsr.bam",
                        help="Preprocessed BAM suffix (default: _bqsr.bam). Use "
                             "_coordsort.bam if preprocess ran without --known-sites")
    parser.add_argument("--germline-resource",
                        help="af-only germline VCF. Optional. ")
    parser.add_argument("--intervals", metavar="FILE",
                        help="GATK -L intervals. Default: every contig in the reference. ")
    parser.add_argument("--merge-contigs-into-num-partitions", type=int, default=25,
                        metavar="N",
                        help="GenomicsDBImport partitioning (default: 25). ")
    parser.add_argument("--max-mnp-distance", type=int, default=0,
                        help="Mutect2 --max-mnp-distance (default: 0). Must be 0 "
                             "for PON construction. ")
    parser.add_argument("--steps", type=parse_steps, default=list(ALL_STEPS),
                        help=f"Comma-separated subset of {','.join(ALL_STEPS)} "
                             f"(default: all three)")
    parser.add_argument("--temp-dir", metavar="DIR", help="GATK temp directory")
    parser.add_argument("--java-heap", default=DEFAULT_JAVA_HEAP, metavar="SIZE",
                        help=f"GATK --java-options -Xmx (default: {DEFAULT_JAVA_HEAP})")
    parser.add_argument("--work-dir", help="Where to write the decompressed FASTA if --ref is .gz")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Control samples to call concurrently in one process (default: 1)")
    parser.add_argument("--sample-index", type=int, metavar="N",
                        help="Call only control sample N (0-based among controls). Set "
                             "automatically from $SLURM_ARRAY_TASK_ID by --array")
    parser.add_argument("--overwrite", action="store_true",
                        help="Recall samples that already have a VCF, and discard an existing "
                             "GenomicsDB workspace")
    parser.add_argument("--hpc", action="store_true", help="Submit as SLURM jobs")
    parser.add_argument("--array", action="store_true",
                        help="With --hpc, run the `call` stage as a job array of one task per "
                             "control, then submit `import,create` as a dependent job")
    parser.add_argument("--array-concurrency", type=int, metavar="N",
                        help="Cap simultaneously running array tasks (SLURM's %%N suffix)")
    parser.add_argument("--conda-env", default="bio_qc",
                        help="Conda environment activated in the generated SLURM script")
    parser.add_argument("--modules", default="", help="Comma-separated HPC modules to load")
    parser.add_argument("--directives-dict", dest="directives_dict",
                        type=parse_directives_dict, default={}, metavar="DICT",
                        help="SLURM sbatch directives as a Python dict literal. Example: "
                             "\"{'partition': 'Mamba', 'mem': '48gb', 'time': '24:00:00'}\"")


def resolve_control_samples(meta, control_group):
    """
    The control samples, with the group name validated against the metadata.

    """
    if control_group not in meta.groups:
        raise SystemExit(
                f"wgseqpipe pon: --control-group {control_group!r} is not one of the "
                f"metadata's group values {meta.groups}. Pick whichever of those two "
                f"is the normal/control condition."
                )
    controls = meta.samples_in_group(control_group)
    if not controls:
        raise SystemExit(f"wgseqpipe pon: no samples in group {control_group!r}")
    return controls


def ensure_intervals(ref_fasta, intervals, output_dir):
    """
    An intervals file, synthesised from the .fai if none was given.

    """
    if intervals:
        return Path(intervals)
    fai = Path(f"{ref_fasta}.fai")
    if not fai.is_file():
        raise SystemExit(
                f"wgseqpipe pon: {fai} not found, so --intervals cannot be derived. Run "
                f"`tarva prep-reference --ref {ref_fasta}` or pass --intervals."
                )
    out = Path(output_dir) / "all_contigs.list"
    if out.is_file():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(fai) as src:
        contigs = [line.split("\t")[0] for line in src if line.strip()]
    out.write_text("\n".join(contigs) + "\n")
    print(f"Wrote {out} ({len(contigs)} contigs from {fai.name})")
    return out


def adopt_existing_pon(existing, output_dir):
    """
    Link an externally built panel into <output>/pon.vcf.gz.

    """
    existing = Path(existing)
    if not existing.is_file():
        raise SystemExit(f"wgseqpipe pon: --existing-pon not found: {existing}")
    index = next((p for p in (Path(f"{existing}.tbi"), Path(f"{existing}.idx"))
                  if p.is_file()), None)
    if index is None:
        raise SystemExit(
                f"wgseqpipe pon: --existing-pon {existing} has no index beside it "
                f"(expected {existing.name}.tbi or .idx). Index it with "
                f"`gatk IndexFeatureFile -I {existing}`."
                )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pon_vcf = output_dir / PON_NAME
    for source, target in ((existing, pon_vcf),
                           (index, output_dir / f"{PON_NAME}{index.suffix}")):
        if target.is_symlink() or target.is_file():
            target.unlink()
        try:
            target.symlink_to(source.resolve())
        except OSError:
            shutil.copy2(source, target)
    (output_dir / "pon_source.txt").write_text(f"{existing.resolve()}\n")
    print(f"Adopted existing panel of normals: {existing}")
    print(f"  linked to {pon_vcf} (provenance recorded in pon_source.txt)")
    print("Pass that path to `wgseqpipe variants --panel-of-normals`.")
    return pon_vcf


def _call_one(sample, args, ref, calls_dir):
    """Mutect2 on one control."""
    out_vcf = calls_dir / f"{sample.individual_id}{PON_VCF_SUFFIX}"
    row = {"individualID": sample.individual_id, "group": sample.group,
           "status": "FAILED", "vcf": "", "error": ""}
    try:
        if out_vcf.is_file() and not args.overwrite:
            print(f"{sample.individual_id}: {out_vcf.name} already present, skipping")
            return {**row, "status": "skipped", "vcf": str(out_vcf)}
        bam = (Path(args.preprocessed) / sample.group / sample.individual_id
               / f"{sample.individual_id}{args.bam_suffix}")
        if not bam.is_file():
            raise FileNotFoundError(
                    f"preprocessed BAM not found: {bam}. Check --preprocessed and "
                    f"--bam-suffix (use _coordsort.bam if preprocess ran without "
                    f"--known-sites)."
                    )
        gatk_args = ["-R", ref, "-I", bam, "-O", out_vcf,
                     "--max-mnp-distance", args.max_mnp_distance]
        if args.germline_resource:
            gatk_args += ["--germline-resource", args.germline_resource]
        run_gatk("Mutect2", gatk_args, temp_dir=args.temp_dir, java_heap=args.java_heap)
        return {**row, "status": "called", "vcf": str(out_vcf)}
    except GatkStepError as exc:
        return {**row, "error": str(exc)}
    except Exception as exc:
        return {**row, "error": f"{type(exc).__name__}: {exc}"}


def _select_controls(controls, sample_index):
    if sample_index is None:
        return controls
    if not 0 <= sample_index < len(controls):
        raise SystemExit(
                f"wgseqpipe pon: --sample-index {sample_index} is out of range for "
                f"{len(controls)} control sample(s) (valid: 0..{len(controls) - 1})."
                )
    return [controls[sample_index]]


def _summary_path(output_dir, controls, sample_index):
    if sample_index is None:
        return Path(output_dir) / "pon_summary.csv"
    return (Path(output_dir) / "summary_parts"
            / f"{sample_index:03d}_{controls[0].individual_id}.csv")


def run_call(args, controls, ref):
    """Stage 1: Mutect2 across the control group."""
    output_dir = Path(args.output)
    calls_dir = output_dir / "calls"
    calls_dir.mkdir(parents=True, exist_ok=True)
    selected = _select_controls(controls, args.sample_index)
    summary_csv = _summary_path(output_dir, selected, args.sample_index)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"PON call: {len(selected)} of {len(controls)} control sample(s) "
          f"({', '.join(s.individual_id for s in selected)}), {args.jobs} at a time")
    rows = []
    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        fh.flush()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(_call_one, s, args, ref, calls_dir)
                       for s in selected]
            for future in as_completed(futures):
                row = future.result()
                writer.writerow(row)
                fh.flush()
                rows.append(row)

    failed = [r for r in rows if r["status"] == "FAILED"]
    print(f"\nPON call: {len(rows) - len(failed)} done, {len(failed)} failed")
    print(f"Wrote summary: {summary_csv}")
    for row in failed:
        print(f"  FAILED {row['individualID']}: {row['error']}", file=sys.stderr)
    if failed:
        raise SystemExit(f"wgseqpipe pon: {len(failed)} control(s) failed -- see {summary_csv}")
    return rows


def run_import(args, controls, ref, intervals):
    """Stage 2: GenomicsDBImport of every control VCF."""
    output_dir = Path(args.output)
    db_path = output_dir / PON_DB_NAME
    vcfs = [output_dir / "calls" / f"{s.individual_id}{PON_VCF_SUFFIX}" for s in controls]
    missing = [str(v) for v in vcfs if not v.is_file()]
    if missing:
        raise SystemExit(
                f"wgseqpipe pon: cannot import, {len(missing)} control VCF(s) missing: "
                f"{missing[:5]}{'...' if len(missing) > 5 else ''}. Run --steps call first."
                )
    if db_path.exists():
        if not args.overwrite:
            raise SystemExit(
                    f"wgseqpipe pon: {db_path} already exists. GenomicsDBImport will not "
                    f"write into an existing workspace -- pass --overwrite to discard it."
                    )
        print(f"--overwrite: removing {db_path}")
        shutil.rmtree(db_path)

    gatk_args = ["-R", ref, "--genomicsdb-workspace-path", db_path, "-L", intervals,
                 "--merge-contigs-into-num-partitions",
                 args.merge_contigs_into_num_partitions]
    for vcf in vcfs:
        gatk_args += ["-V", vcf]
    run_gatk("GenomicsDBImport", gatk_args, temp_dir=args.temp_dir, java_heap=args.java_heap)
    print(f"Wrote GenomicsDB workspace: {db_path}")
    return db_path


def run_create(args, ref):
    """
    Stage 3: CreateSomaticPanelOfNormals, straight from the GenomicsDB.

    """
    output_dir = Path(args.output)
    db_path = output_dir / PON_DB_NAME
    if not db_path.exists():
        raise SystemExit(
                f"wgseqpipe pon: {db_path} not found. Run --steps import first."
                )
    pon_vcf = output_dir / PON_NAME
    # No --min-sample-count: legacy passed none and took GATK's default of 2,
    # so emitting nothing reproduces it exactly.
    gatk_args = ["-R", ref, "-V", f"gendb://{db_path}", "-O", pon_vcf]
    if args.germline_resource:
        gatk_args += ["--germline-resource", args.germline_resource]
    run_gatk("CreateSomaticPanelOfNormals", gatk_args,
             temp_dir=args.temp_dir, java_heap=args.java_heap)
    print(f"\nWrote panel of normals: {pon_vcf}")
    print("Pass this to `wgseqpipe variants --panel-of-normals`, which applies it to "
          "BOTH groups.")
    return pon_vcf


def run(args):
    if args.existing_pon:
        adopt_existing_pon(args.existing_pon, args.output)
        return
    if args.hpc:
        _submit_hpc(args)
        return
    meta = load_metadata(args.metadata)
    controls = resolve_control_samples(meta, args.control_group)
    other = [g for g in meta.groups if g != args.control_group]
    print(f"Control group: {args.control_group!r} ({len(controls)} samples); "
          f"other group: {other[0]!r}")
    ref = ensure_unzipped_fasta(args.ref, work_dir=args.work_dir)
    intervals = ensure_intervals(ref, args.intervals, args.output)
    if not args.germline_resource:
        print("No --germline-resource: Mutect2 will not suppress common sites. ", file=sys.stderr)

    if "call" in args.steps:
        run_call(args, controls, ref)
    if "import" in args.steps:
        run_import(args, controls, ref, intervals)
    if "create" in args.steps:
        run_create(args, ref)


def _base_command(args, steps):
    parts = ["python", "-u", "-m", "tarva.cli", "wgseqpipe", "pon",
             "--metadata", str(args.metadata),
             "--preprocessed", str(args.preprocessed),
             "--output", str(args.output),
             "--ref", str(args.ref),
             "--control-group", str(args.control_group),
             "--bam-suffix", str(args.bam_suffix),
             "--merge-contigs-into-num-partitions",
             str(args.merge_contigs_into_num_partitions),
             "--max-mnp-distance", str(args.max_mnp_distance),
             "--java-heap", str(args.java_heap),
             "--jobs", str(args.jobs),
             "--steps", ",".join(steps)]
    for flag, value in (("--germline-resource", args.germline_resource),
                        ("--intervals", args.intervals),
                        ("--temp-dir", args.temp_dir),
                        ("--work-dir", args.work_dir)):
        if value:
            parts += [flag, str(value)]
    if args.overwrite:
        parts.append("--overwrite")
    return parts


def _write_and_submit(args, steps, job_name, sbatch_kwargs, extra_directives,
                      command_parts, log_suffix=""):
    log_dir = Path(args.output) / "logs"
    script_path = write_sbatch_script(
            job_name=job_name,
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / f"{job_name}.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            log_suffix=log_suffix,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    output = submit_job(script_path)
    print(output)
    match = re.search(r"(\d+)", output)
    return match.group(1) if match else None


def _submit_hpc(args):
    """
    Submit the requested stages.

    """
    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    sbatch_kwargs.setdefault("cpus_per_task", 4 * args.jobs)

    if not args.array:
        _write_and_submit(args, args.steps, "tarva_wgs_pon", sbatch_kwargs,
                          extra_directives, _base_command(args, args.steps))
        return

    meta = load_metadata(args.metadata)
    controls = resolve_control_samples(meta, args.control_group)
    array_job = None

    if "call" in args.steps:
        array_spec = f"0-{len(controls) - 1}"
        if args.array_concurrency:
            array_spec += f"%{args.array_concurrency}"
        call_directives = dict(extra_directives)
        call_directives.setdefault("array", array_spec)
        command = _base_command(args, ["call"]) + ["--sample-index", "$SLURM_ARRAY_TASK_ID"]
        print(f"Submitting PON call array of {len(controls)} task(s): --array={array_spec}")
        array_job = _write_and_submit(args, ["call"], "tarva_wgs_pon_call", sbatch_kwargs,
                                      call_directives, command, log_suffix="_%a")

    aggregate = [s for s in args.steps if s in ("import", "create")]
    if not aggregate:
        return
    agg_directives = dict(extra_directives)
    if array_job:
        # Appended rather than setdefault
        existing = agg_directives.get("dependency")
        agg_directives["dependency"] = (f"{existing},afterok:{array_job}" if existing
                                        else f"afterok:{array_job}")
        print(f"Submitting import/create with dependency "
              f"{agg_directives['dependency']}")
    # One process
    agg_kwargs = dict(sbatch_kwargs, cpus_per_task=4)
    _write_and_submit(args, aggregate, "tarva_wgs_pon_build", agg_kwargs,
                      agg_directives, _base_command(args, aggregate))
