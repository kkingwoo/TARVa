"""
analysispipe validate: independent checks on the estimated proportions.

    `tarva analysispipe validate --output <results_dir> [--label 5type_tmm]
        [--expression tmm] [--aggregator median] [--n-marker 100]
        [--hpc] [--r-module R/4.5.1] [--r-libs PATH] [--directives-dict "{...}"]`

    Input:  the matrices and proportions `expression` and `dsa` wrote into the
            same --output.
    Output: 07_marker_check.csv, 07_svd_comparison.csv,
            07_pc1_comparison.csv, 07_validate_summary.json, 07_validation.png

    Deconvolution fails quietly, so no single method's output should be
    presented unchecked. Three checks, each able to fail independently. All
    three are INTERNAL -- they compare this run against other quantities
    derived from the same cohort, so a disagreement localises the problem
    rather than merely flagging one:

      1. Direct marker expression by group. If SNAP25 is flat while DSA reports
         a large neuronal shift, something upstream is wrong -- this needs no
         algorithm to interpret.

      2. BRETIGEA's own SVD estimator on the same markers. Different
         mathematics, so agreement is real corroboration. Its output is
         surrogate variables on an arbitrary scale, so compare Spearman rank
         correlation across samples, never absolute values.

      3. PC1 of the expression matrix against the neuronal proportion. In bulk
         cortex the dominant axis of variation is usually neuronal content; if
         the two track, two independent routes agree on what dominates the
         transcriptome.

    There is deliberately no check against literature proportions. A
    fold-change against a number from a different cohort, region and method
    reads as a validation while testing nothing; such comparisons belong in
    prose, where the provenance of each figure can be stated.
"""

from pathlib import Path

from tarva.analysispipe.rstep import add_common_arguments, run_r_step

DEFAULT_CHECK_GENES = (
    "SNAP25,RBFOX3,SYT1,GFAP,AQP4,SLC1A2,CD68,AIF1,CSF1R,P2RY12,"
    "MBP,PLP1,MOBP,CLDN5,FLT1,PECAM1"
)


def _csv_list(raw):
    return [v.strip() for v in str(raw).split(",") if v.strip()]


def add_arguments(parser):
    parser.add_argument(
            "--label", default=None, metavar="NAME",
            help="Label of the dsa run to validate, matching that run's --label",
            )
    parser.add_argument(
            "--expression", default=None, choices=("tpm", "tmm"),
            help=(
                "Expression matrix to validate against. Defaults to whatever the dsa "
                "run used, read from its summary."
                ),
            )
    parser.add_argument(
            "--aggregator", default="median",
            help="Aggregator to validate (default: median)",
            )
    parser.add_argument(
            "--check-genes", default=DEFAULT_CHECK_GENES,
            help="Comma-separated canonical markers for the direct expression check",
            )
    parser.add_argument(
            "--n-marker", type=int, default=100,
            help="Markers per cell type for BRETIGEA's SVD estimator (default: 100)",
            )
    add_common_arguments(parser)


def run(args):
    output_dir = Path(args.output)
    sfx = f"_{args.label}" if args.label else ""
    for needed in ("02_expression.rds", f"05_cell_proportions{sfx}.csv"):
        if not (output_dir / needed).is_file():
            raise SystemExit(
                    f"tarva analysispipe: {output_dir / needed} not found. Run the "
                    f"expression and dsa steps into {output_dir} first."
                    )
    run_r_step(
            args,
            step="validate",
            script_name="07_validate.R",
            extra={
                "label":       args.label,
                "expression":  args.expression,
                "aggregator":  args.aggregator,
                "check_genes": _csv_list(args.check_genes),
                "n_marker":    args.n_marker,
                },
            job_name=f"tarva_validate{sfx}",
            )
