"""
analysispipe coherence: marker co-expression coherence per cell type.

    `tarva analysispipe coherence --output <results_dir>
        [--top-n 500] [--coherence-min 0.25]
        [--hpc] [--r-module R/4.5.1] [--r-libs PATH] [--directives-dict "{...}"]`

    Input:  the expression matrices `expression` wrote into the same --output.
    Output: 04_marker_coherence.csv, 04_coherence_summary.json, 04_coherence.png

    A cell type's markers should rise and fall together across samples, driven
    by that cell type's abundance. Median pairwise Pearson correlation among a
    type's markers measures whether they do. DSA assumes marker average
    expression is proportional to cell-type proportion; if the markers do not
    co-vary, that average tracks nothing in particular, and because nnls solves
    for all cell types jointly an incoherent type can also distort the others.

    Every cell type in the marker table is measured, so any exclusion is
    documented from this dataset rather than inherited. Types falling below
    --coherence-min are reported as failing, and `dsa` uses the passing set by
    default.

    Correlations are computed on log2 because expression has a heavy right tail
    that would otherwise let a few high-expressing genes dominate Pearson. DSA
    itself still uses the linear TPM.
"""

from pathlib import Path

from tarva.analysispipe.rstep import add_common_arguments, run_r_step


def add_arguments(parser):
    parser.add_argument(
            "--top-n", type=int, default=500,
            help="Rank-ordered markers per cell type to use (default: 500)",
            )
    parser.add_argument(
            "--coherence-min", type=float, default=0.25,
            help=(
                "Median pairwise correlation below which a cell type is reported as "
                "failing and excluded from step 5 by default (default: 0.25). Wang et "
                "al. found OPC at 0.06-0.12 and the retained five at 0.4-0.8."
                ),
            )
    add_common_arguments(parser)


def run(args):
    output_dir = Path(args.output)
    if not (output_dir / "02_expression.rds").is_file():
        raise SystemExit(
                f"tarva analysispipe: {output_dir / '02_expression.rds'} not found. "
                f"Run `tarva analysispipe expression --output {output_dir}` first."
                )
    run_r_step(
            args,
            step="coherence",
            script_name="04_marker_coherence.R",
            extra={"top_n": args.top_n, "coherence_min": args.coherence_min},
            )
