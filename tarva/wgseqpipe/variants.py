"""
wgseqpipe variants: call and filter WGS variants for every sample, applying
the panel of normals to both groups.

Per sample: Mutect2 with the panel of normals, then GetPileupSummaries ->
CalculateContamination -> FilterMutectCalls.

    `tarva wgseqpipe variants --metadata <meta.csv> --preprocessed <preproc_dir>
        --output <vcf_dir> --ref <genome.fa[.gz]>
        [--panel-of-normals <pon.vcf.gz>] [--germline-resource <af-only.vcf.gz>]
        [--common-sites <small_exac_common.vcf.gz>]
        [--max-mnp-distance N] [--java-heap SIZE] [--temp-dir DIR]
        [--jobs N] [--sample-index N] [--overwrite]
        [--hpc] [--array] [--array-concurrency N] [--directives-dict "{...}"]`

Input:  the preprocessed BAMs and the reference; optionally a panel of
        normals, a germline resource and a common-sites VCF. Without
        --common-sites the pileup and contamination steps are skipped and the
        sample is reported as filtered with no contamination table.
Output: <output>/<group>/<individualID>/<individualID>_contaminates_filtered.vcf.gz
        per sample, the raw Mutect2 VCF beside it, and a summary CSV. This is
        what `analysispipe rna-unique` reads as the WGS side.
"""

import csv
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tarva.gatk_common import DEFAULT_JAVA_HEAP, GatkStepError, run_gatk
from tarva.metadata import load_metadata
from tarva.reference import ensure_unzipped_fasta
from tarva.slurm_utils import (parse_directives_dict, split_directives, submit_job,
                               write_sbatch_script)

RAW_VCF_SUFFIX = "_mutect2.vcf.gz"
# The name `analysispipe rna-unique` looks for. Keep it.
FILTERED_VCF_SUFFIX = "_contaminates_filtered.vcf.gz"

SUMMARY_FIELDS = [
    "individualID", "group", "status", "filtered_vcf", "raw_vcf",
    "contamination_table", "contamination", "error",
    ]


def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--preprocessed", required=True,
                        help="The --output directory `wgseqpipe preprocess` wrote")
    parser.add_argument("--output", required=True,
                        help="Directory for VCFs, as <output>/<group>/<individualID>/")
    parser.add_argument("--ref", required=True,
                        help="Reference genome FASTA; must match the preprocessed BAMs")
    parser.add_argument("--bam-suffix", default="_bqsr.bam",
                        help="Preprocessed BAM suffix (default: _bqsr.bam). Use "
                             "_coordsort.bam if preprocess ran without --known-sites")
    parser.add_argument("--panel-of-normals", metavar="VCF",
                        help="pon.vcf.gz from `wgseqpipe pon`. Optional. ")
    parser.add_argument("--germline-resource", metavar="VCF",
                        help="af-only germline VCF. Optional. ")
    parser.add_argument("--common-sites", metavar="VCF",
                        help="Common biallelic sites for GetPileupSummaries. Optional. ")
    parser.add_argument("--max-mnp-distance", type=int, default=0,
                        help="Mutect2 --max-mnp-distance (default: 0) ")
    parser.add_argument("--temp-dir", metavar="DIR", help="GATK temp directory")
    parser.add_argument("--java-heap", default=DEFAULT_JAVA_HEAP, metavar="SIZE",
                        help=f"GATK --java-options -Xmx (default: {DEFAULT_JAVA_HEAP})")
    parser.add_argument("--work-dir", help="Where to write the decompressed FASTA if --ref is .gz")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Samples to call concurrently in one process (default: 1)")
    parser.add_argument("--sample-index", type=int, metavar="N",
                        help="Call only sample N (0-based, metadata row order). Set "
                             "automatically from $SLURM_ARRAY_TASK_ID by --array")
    parser.add_argument("--overwrite", action="store_true",
                        help="Recall samples that already have a filtered VCF")
    parser.add_argument("--hpc", action="store_true", help="Submit as a SLURM job")
    parser.add_argument("--array", action="store_true",
                        help="With --hpc, submit a job array of one task per sample")
    parser.add_argument("--array-concurrency", type=int, metavar="N",
                        help="Cap simultaneously running array tasks (SLURM's %%N suffix)")
    parser.add_argument("--conda-env", default="bio_qc",
                        help="Conda environment activated in the generated SLURM script")
    parser.add_argument("--modules", default="", help="Comma-separated HPC modules to load")
    parser.add_argument("--directives-dict", dest="directives_dict",
                        type=parse_directives_dict, default={}, metavar="DICT",
                        help="SLURM sbatch directives as a Python dict literal. Example: "
                             "\"{'partition': 'Mamba', 'mem': '48gb', 'time': '24:00:00'}\"")


def _contamination_fraction(table):
    """The fraction from CalculateContamination's two-line table, for the summary."""
    try:
        with open(table) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return rows[0].get("contamination", "") if rows else ""
    except (OSError, csv.Error, IndexError):
        return ""


def _variants_one(sample, args, ref):
    """
    One sample at a time.
    """
    sample_dir = Path(args.output) / sample.group / sample.individual_id
    raw_vcf = sample_dir / f"{sample.individual_id}{RAW_VCF_SUFFIX}"
    filtered_vcf = sample_dir / f"{sample.individual_id}{FILTERED_VCF_SUFFIX}"
    pileup_table = sample_dir / f"{sample.individual_id}_pileupSum.table"
    contam_table = sample_dir / f"{sample.individual_id}_contamination.table"
    row = {"individualID": sample.individual_id, "group": sample.group,
           "status": "FAILED", "filtered_vcf": "", "raw_vcf": "",
           "contamination_table": "", "contamination": "", "error": ""}
    try:
        if filtered_vcf.is_file() and not args.overwrite:
            print(f"{sample.individual_id}: {filtered_vcf.name} already present, skipping")
            return {**row, "status": "skipped", "filtered_vcf": str(filtered_vcf),
                    "raw_vcf": str(raw_vcf) if raw_vcf.is_file() else ""}

        bam = (Path(args.preprocessed) / sample.group / sample.individual_id
               / f"{sample.individual_id}{args.bam_suffix}")
        if not bam.is_file():
            raise FileNotFoundError(
                    f"preprocessed BAM not found: {bam}. Check --preprocessed and "
                    f"--bam-suffix."
                    )
        sample_dir.mkdir(parents=True, exist_ok=True)

        mutect_args = ["-R", ref, "-I", bam, "-O", raw_vcf,
                       "--max-mnp-distance", args.max_mnp_distance]
        if args.panel_of_normals:
            mutect_args += ["--panel-of-normals", args.panel_of_normals]
        if args.germline_resource:
            mutect_args += ["--germline-resource", args.germline_resource]
        # Resume point: Mutect2 is the expensive step, so a failure in the
        # contamination or filtering stages below should not discard it.
        if raw_vcf.is_file() and not args.overwrite:
            print(f"{sample.individual_id}: {raw_vcf.name} already present, skipping Mutect2")
        else:
            run_gatk("Mutect2", mutect_args, temp_dir=args.temp_dir, java_heap=args.java_heap)
        row["raw_vcf"] = str(raw_vcf)

        filter_args = ["-R", ref, "-V", raw_vcf, "-O", filtered_vcf]
        status = "filtered"
        if args.common_sites:
            # -V and -L are both the common-sites VCF.
            run_gatk("GetPileupSummaries",
                     ["-I", bam, "-V", args.common_sites, "-L", args.common_sites,
                      "-O", pileup_table],
                     temp_dir=args.temp_dir, java_heap=args.java_heap)
            run_gatk("CalculateContamination",
                     ["-I", pileup_table, "-O", contam_table],
                     temp_dir=args.temp_dir, java_heap=args.java_heap)
            filter_args += ["--contamination-table", contam_table]
            row["contamination_table"] = str(contam_table)
            row["contamination"] = _contamination_fraction(contam_table)
        else:
            status = "filtered (no contamination)"
        run_gatk("FilterMutectCalls", filter_args,
                 temp_dir=args.temp_dir, java_heap=args.java_heap)
        return {**row, "status": status, "filtered_vcf": str(filtered_vcf)}
    except GatkStepError as exc:
        return {**row, "error": str(exc)}
    except Exception as exc:
        return {**row, "error": f"{type(exc).__name__}: {exc}"}


def _select_samples(meta, sample_index):
    if sample_index is None:
        return meta.samples
    if not 0 <= sample_index < len(meta.samples):
        raise SystemExit(
                f"wgseqpipe variants: --sample-index {sample_index} is out of range for "
                f"{len(meta.samples)} sample(s) (valid: 0..{len(meta.samples) - 1})."
                )
    return [meta.samples[sample_index]]


def _summary_path(output_dir, samples, sample_index):
    if sample_index is None:
        return Path(output_dir) / "wgs_variants_summary.csv"
    return (Path(output_dir) / "summary_parts"
            / f"{sample_index:03d}_{samples[0].individual_id}.csv")


def _merge_summary_parts(output_dir, total_samples):
    parts = sorted((Path(output_dir) / "summary_parts").glob("*.csv"))
    combined = Path(output_dir) / "wgs_variants_summary.csv"
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
                for sub_row in csv.DictReader(part_fh):
                    writer.writerow(sub_row)
    partial.replace(combined)
    print(f"Wrote combined summary: {combined}")
    return combined


def run(args):
    if args.hpc:
        _submit_hpc(args)
        return
    meta = load_metadata(args.metadata)
    samples = _select_samples(meta, args.sample_index)
    ref = ensure_unzipped_fasta(args.ref, work_dir=args.work_dir)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = _summary_path(output_dir, samples, args.sample_index)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reference: {ref}")
    print(f"Mutect2: {len(samples)} of {len(meta.samples)} sample(s) "
          f"({', '.join(s.individual_id for s in samples)}), {args.jobs} at a time")
    print(f"Panel of normals: {args.panel_of_normals or 'none'} -- applied to "
          f"{'all groups' if args.panel_of_normals else 'nothing'}")
    print(f"Germline resource: {args.germline_resource or 'none'}")
    if not args.common_sites:
        print("No --common-sites: GetPileupSummaries and CalculateContamination are "
              "skipped and FilterMutectCalls runs without a contamination table.",
              file=sys.stderr)

    rows = []
    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        fh.flush()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(_variants_one, s, args, ref) for s in samples]
            for future in as_completed(futures):
                row = future.result()
                writer.writerow(row)
                fh.flush()
                rows.append(row)

    failed = [r for r in rows if r["status"] == "FAILED"]
    print(f"\nWGS variants: {len(rows) - len(failed)} done, {len(failed)} failed")
    print(f"Wrote summary: {summary_csv}")
    if args.sample_index is not None:
        _merge_summary_parts(output_dir, len(meta.samples))
    for row in failed:
        print(f"  FAILED {row['individualID']}: {row['error']}", file=sys.stderr)
    if failed:
        raise SystemExit(
                f"wgseqpipe variants: {len(failed)} sample(s) failed -- see {summary_csv}"
                )
    return rows


def _submit_hpc(args):
    command_parts = ["python", "-u", "-m", "tarva.cli", "wgseqpipe", "variants",
                     "--metadata", str(args.metadata),
                     "--preprocessed", str(args.preprocessed),
                     "--output", str(args.output),
                     "--ref", str(args.ref),
                     "--bam-suffix", str(args.bam_suffix),
                     "--max-mnp-distance", str(args.max_mnp_distance),
                     "--java-heap", str(args.java_heap),
                     "--jobs", str(args.jobs)]
    for flag, value in (("--panel-of-normals", args.panel_of_normals),
                        ("--germline-resource", args.germline_resource),
                        ("--common-sites", args.common_sites),
                        ("--temp-dir", args.temp_dir),
                        ("--work-dir", args.work_dir)):
        if value:
            command_parts += [flag, str(value)]
    if args.overwrite:
        command_parts.append("--overwrite")

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    sbatch_kwargs.setdefault("cpus_per_task", 4 * args.jobs)

    log_suffix = ""
    if args.array:
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
            job_name="tarva_wgs_variants",
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / "tarva_wgs_variants.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            log_suffix=log_suffix,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
