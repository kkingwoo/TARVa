"""
rnaseqpipe variants: HaplotypeCaller and filtering for RNA variant calls

    `tarva rnaseqpipe variants --metadata <meta.csv> --preprocessed <preproc_dir>
        --output <vcf_dir> --ref <reference.fa[.gz]>
        [--bam-suffix SUFFIX] [--output-mode MODE] [--analyze-covariates]
        [--qual-filter F] [--qual-depth-filter F] [--depth-filter N]
        [--java-heap SIZE] [--temp-dir DIR] [--jobs N] [--sample-index N]
        [--overwrite] [--hpc] [--array] [--array-concurrency N]
        [--directives-dict "{...}"]`

Per sample: optional AnalyzeCovariates, then HaplotypeCaller, then
VariantFiltration. Output is
<output>/<group>/<individualID>/<individualID>_varcall_flt.vcf, which is the
name TARVaCreation/OriginalBuild/checkRNA_againstWGS.py tests for with
`f.endswith("varcall_flt.vcf")`.

**That name is uncompressed on purpose.** The endswith test does not match
`.vcf.gz`, so bgzipping this output silently makes the RNA half of the WGS/RNA
comparison invisible. Legacy named it `_filtered.vcf` and so never matched
either -- a pre-existing broken handoff this module fixes by naming the file
what the consumer actually looks for.

Replaces the calling half of legacy GatkRnaSeqPipe/gatkPipe.py (haplocall,
filter_vars, analyzecovariates). The filter expressions are legacy's, unchanged
and exposed as flags rather than hardcoded:

    QUAL < 25.0   -> qual_filt
    QD   < 10.0   -> qual_depth_filt
    DP   < 20     -> depth_filt

The HaplotypeCaller and VariantFiltration argument lists are legacy's exactly.
No -L, and none of GATK's RNA-specific options.

One deviation: `--analyze-covariates` is off by default where legacy called
AnalyzeCovariates unconditionally. It renders its plots through R, and section 2
of the status doc records that loading the R module alongside an active conda
env breaks R via LD_LIBRARY_PATH shadowing. It produces a diagnostic PDF, not
pipeline output, so a failure there should not take the run down with it. When
enabled, a failure is recorded in the summary and the sample continues.

`--output-mode` defaults to EMIT_ALL_CONFIDENT_SITES, as legacy. Be aware that
this emits confidently-called reference sites too, not just variants, so the VCF
is far larger than a variants-only call -- which is the point, since the
downstream comparison asks whether a position was callable at all, not only
whether it varied.
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

RAW_VCF_SUFFIX = "_haplo.vcf"
# checkRNA_againstWGS.py:31 does f.endswith("varcall_flt.vcf") -- not .gz.
FILTERED_VCF_SUFFIX = "_varcall_flt.vcf"
DEFAULT_OUTPUT_MODE = "EMIT_ALL_CONFIDENT_SITES"
DEFAULT_QUAL_FILTER = 25.0
DEFAULT_QUAL_DEPTH_FILTER = 10.0
DEFAULT_DEPTH_FILTER = 20

SUMMARY_FIELDS = [
    "individualID", "group", "status", "filtered_vcf", "raw_vcf",
    "covariates_pdf", "warning", "error",
    ]


def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--preprocessed", required=True,
                        help="The --output directory `rnaseqpipe preprocess` wrote")
    parser.add_argument("--output", required=True,
                        help="Directory for VCFs, as <output>/<group>/<individualID>/")
    parser.add_argument("--ref", required=True,
                        help="Reference FASTA; must match the preprocessed BAMs -- normally "
                             "the combined_reference.fa the two-pass alignment used")
    parser.add_argument("--bam-suffix", default="_bqsr.bam",
                        help="Preprocessed BAM suffix (default: _bqsr.bam). Use _rg.bam if "
                             "preprocess ran without --known-sites, or _splitncigar.bam if "
                             "it also ran with --read-group-repair never")
    parser.add_argument("--output-mode", default=DEFAULT_OUTPUT_MODE,
                        help=f"HaplotypeCaller --output-mode (default: {DEFAULT_OUTPUT_MODE}, "
                             "as legacy). EMIT_VARIANTS_ONLY gives a far smaller VCF but "
                             "loses the ability to distinguish 'reference here' from 'not "
                             "callable here' in the WGS/RNA comparison")
    parser.add_argument("--analyze-covariates", action="store_true",
                        help="Also run AnalyzeCovariates on the BQSR table from preprocess, "
                             "producing a diagnostic PDF. Off by default: it renders through "
                             "R, and an R module loaded beside a conda env breaks on "
                             "LD_LIBRARY_PATH shadowing. Failures here are reported, not fatal")
    parser.add_argument("--qual-filter", type=float, default=DEFAULT_QUAL_FILTER,
                        help=f"VariantFiltration QUAL threshold (default: {DEFAULT_QUAL_FILTER}, "
                             f"as legacy); records below it get FILTER=qual_filt")
    parser.add_argument("--qual-depth-filter", type=float, default=DEFAULT_QUAL_DEPTH_FILTER,
                        help=f"VariantFiltration QD threshold (default: "
                             f"{DEFAULT_QUAL_DEPTH_FILTER}, as legacy) -> qual_depth_filt")
    parser.add_argument("--depth-filter", type=int, default=DEFAULT_DEPTH_FILTER,
                        help=f"VariantFiltration DP threshold (default: {DEFAULT_DEPTH_FILTER}, "
                             f"as legacy) -> depth_filt")
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
                             "\"{'partition': 'Orion', 'mem': '48gb', 'time': '24:00:00'}\"")


def filtration_args(args):
    """
    Legacy's three filters, as (name, expression) pairs.

    VariantFiltration pairs --filter-name with --filter-expression positionally,
    so the two lists must stay in lockstep -- hence building them together here
    rather than appending to one command list in two places.
    """
    pairs = [
        ("qual_filt", f"QUAL < {args.qual_filter}"),
        ("qual_depth_filt", f"QD < {args.qual_depth_filter}"),
        ("depth_filt", f"DP < {args.depth_filter}"),
        ]
    flags = []
    for name, expression in pairs:
        flags += ["--filter-name", name, "--filter-expression", expression]
    return flags


def _variants_one(sample, args, ref):
    sample_dir = Path(args.output) / sample.group / sample.individual_id
    raw_vcf = sample_dir / f"{sample.individual_id}{RAW_VCF_SUFFIX}"
    filtered_vcf = sample_dir / f"{sample.individual_id}{FILTERED_VCF_SUFFIX}"
    row = {"individualID": sample.individual_id, "group": sample.group,
           "status": "FAILED", "filtered_vcf": "", "raw_vcf": "",
           "covariates_pdf": "", "warning": "", "error": ""}
    try:
        if filtered_vcf.is_file() and not args.overwrite:
            print(f"{sample.individual_id}: {filtered_vcf.name} already present, skipping")
            return {**row, "status": "skipped", "filtered_vcf": str(filtered_vcf)}

        sample_bam_dir = Path(args.preprocessed) / sample.group / sample.individual_id
        bam = sample_bam_dir / f"{sample.individual_id}{args.bam_suffix}"
        if not bam.is_file():
            raise FileNotFoundError(
                    f"preprocessed BAM not found: {bam}. Check --preprocessed and "
                    f"--bam-suffix (use _rg.bam if preprocess ran without --known-sites)."
                    )
        sample_dir.mkdir(parents=True, exist_ok=True)

        if args.analyze_covariates:
            recal_table = sample_bam_dir / f"{sample.individual_id}_recal.table"
            pdf = sample_dir / f"{sample.individual_id}_covars_analyzed.pdf"
            if not recal_table.is_file():
                row["warning"] = (f"--analyze-covariates given but {recal_table.name} is "
                                  f"missing (preprocess ran without --known-sites?)")
                print(f"WARNING: {sample.individual_id}: {row['warning']}", file=sys.stderr)
            else:
                # Diagnostic only, and its R dependency is the fragile part of
                # this environment, so a failure must not lose the call set.
                try:
                    run_gatk("AnalyzeCovariates", ["-bqsr", recal_table, "-plots", pdf],
                             temp_dir=args.temp_dir, java_heap=args.java_heap)
                    row["covariates_pdf"] = str(pdf)
                except GatkStepError as exc:
                    row["warning"] = f"AnalyzeCovariates failed ({exc}); continuing"
                    print(f"WARNING: {sample.individual_id}: {row['warning']}",
                          file=sys.stderr)

        # Argument-for-argument what gatkPipe.py's haplocall ran. No -L, and
        # none of GATK's RNA-specific options: legacy passed none of them.
        hc_args = ["-R", ref, "-I", bam, "-O", raw_vcf,
                   "--output-mode", args.output_mode]
        # Resume point: HaplotypeCaller is the expensive step, so a failure in
        # VariantFiltration below should not discard it.
        if raw_vcf.is_file() and not args.overwrite:
            print(f"{sample.individual_id}: {raw_vcf.name} already present, "
                  f"skipping HaplotypeCaller")
        else:
            run_gatk("HaplotypeCaller", hc_args, temp_dir=args.temp_dir,
                     java_heap=args.java_heap)
        row["raw_vcf"] = str(raw_vcf)

        run_gatk("VariantFiltration",
                 ["-R", ref, "-V", raw_vcf, "-O", filtered_vcf] + filtration_args(args),
                 temp_dir=args.temp_dir, java_heap=args.java_heap)
        return {**row, "status": "filtered", "filtered_vcf": str(filtered_vcf)}
    except GatkStepError as exc:
        return {**row, "error": str(exc)}
    except Exception as exc:
        return {**row, "error": f"{type(exc).__name__}: {exc}"}


def _select_samples(meta, sample_index):
    if sample_index is None:
        return meta.samples
    if not 0 <= sample_index < len(meta.samples):
        raise SystemExit(
                f"rnaseqpipe variants: --sample-index {sample_index} is out of range for "
                f"{len(meta.samples)} sample(s) (valid: 0..{len(meta.samples) - 1})."
                )
    return [meta.samples[sample_index]]


def _summary_path(output_dir, samples, sample_index):
    if sample_index is None:
        return Path(output_dir) / "rna_variants_summary.csv"
    return (Path(output_dir) / "summary_parts"
            / f"{sample_index:03d}_{samples[0].individual_id}.csv")


def _merge_summary_parts(output_dir, total_samples):
    parts = sorted((Path(output_dir) / "summary_parts").glob("*.csv"))
    combined = Path(output_dir) / "rna_variants_summary.csv"
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
    print(f"HaplotypeCaller: {len(samples)} of {len(meta.samples)} sample(s) "
          f"({', '.join(s.individual_id for s in samples)}), {args.jobs} at a time")
    print(f"--output-mode {args.output_mode}; filters QUAL<{args.qual_filter}, "
          f"QD<{args.qual_depth_filter}, DP<{args.depth_filter}")

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
    print(f"\nRNA variants: {len(rows) - len(failed)} done, {len(failed)} failed")
    print(f"Wrote summary: {summary_csv}")
    if args.sample_index is not None:
        _merge_summary_parts(output_dir, len(meta.samples))
    for row in failed:
        print(f"  FAILED {row['individualID']}: {row['error']}", file=sys.stderr)
    if failed:
        raise SystemExit(
                f"rnaseqpipe variants: {len(failed)} sample(s) failed -- see {summary_csv}"
                )
    return rows


def _submit_hpc(args):
    command_parts = ["python", "-u", "-m", "tarva.cli", "rnaseqpipe", "variants",
                     "--metadata", str(args.metadata),
                     "--preprocessed", str(args.preprocessed),
                     "--output", str(args.output),
                     "--ref", str(args.ref),
                     "--bam-suffix", str(args.bam_suffix),
                     "--output-mode", str(args.output_mode),
                     "--qual-filter", str(args.qual_filter),
                     "--qual-depth-filter", str(args.qual_depth_filter),
                     "--depth-filter", str(args.depth_filter),
                     "--java-heap", str(args.java_heap),
                     "--jobs", str(args.jobs)]
    for flag, value in (("--temp-dir", args.temp_dir),
                        ("--work-dir", args.work_dir)):
        if value is not None:
            command_parts += [flag, str(value)]
    if args.analyze_covariates:
        command_parts.append("--analyze-covariates")
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
            job_name="tarva_rna_variants",
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / "tarva_rna_variants.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            log_suffix=log_suffix,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
