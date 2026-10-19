"""
analysispipe dsa: Digital Sorting Algorithm cell-type proportions, with
    Pottier's marker-stability procedure.

    `tarva analysispipe dsa --output <results_dir>
        [--cell-types ast,end,mic,neu,oli] [--top-n 500]
        [--n-markers 100] [--n-iter 100] [--seed 20260904]
        [--aggregators median,mean,scaled_mean] [--stability-min 0.8]
        [--hpc] [--r-module R/4.5.1] [--r-libs PATH] [--directives-dict "{...}"]`

    Input:  the matrices and coherence results `expression` and `coherence`
            wrote into the same --output. Cell types default to those passing
            the coherence threshold, so any exclusion follows from the data.
    Output: 05_cell_proportions.csv (long, keyed on rids),
            05_dsa_stability.csv, 05_dsa_markers.csv, 05_dsa_summary.json,
            05_proportions.png. These are RELATIVE proportions of the tissue
            mRNA pool, not cell counts -- neurons carry far more mRNA per cell
            than microglia.

    DSA (Zhong et al. 2013) estimates relative cell-type proportions from
    marker average expression alone, with no reference signature matrix. For
    each cell type it takes the average expression of that type's markers per
    sample, then solves by non-negative least squares for the per-type scaling
    constants that make the weighted sum equal 1 in every sample.

    Median rather than mean aggregation follows Wang et al. 2020: the mean is
    dominated by the few highest-expressed markers, and GFAP in particular is
    upregulated in reactive astrocytes independently of astrocyte number, which
    would confound astrogliosis with astrocyte proliferation. Proportionality
    survives the median because a positive scalar passes through it. Mean and
    per-gene-scaled mean are computed as sensitivity checks, not alternatives.

    Following Pottier et al. 2022, --n-markers are sampled from the top --top-n
    per cell type and proportions estimated across --n-iter iterations;
    between-run Pearson correlations assess stability and the run correlating
    best with all others is retained as representative.

    The reported proportions are the MEAN over those iterations when every
    retained cell type clears --stability-min, with the between-iteration SD
    written alongside them -- a single draw discards the ensemble the stability
    loop paid for, and downstream steps treat these values as known covariates,
    so their sampling error belongs on the record. When a type fails the
    stability gate the estimator falls back to the representative run and says
    so; see --proportions.

"""

from pathlib import Path

from tarva.analysispipe.rstep import add_common_arguments, run_r_step

DEFAULT_AGGREGATORS = "median,mean,scaled_mean"


def _csv_list(raw):
    return [v.strip() for v in str(raw).split(",") if v.strip()]


def add_arguments(parser):
    parser.add_argument(
            "--cell-types", default=None,
            help=(
                "Comma-separated cell types to deconvolve. Default: those passing "
                "step 4's coherence threshold, read from 04_marker_coherence.csv."
                ),
            )
    parser.add_argument(
            "--top-n", type=int, default=500,
            help="Marker pool per cell type to sample from (default: 500)",
            )
    parser.add_argument(
            "--n-markers", type=int, default=100,
            help="Markers sampled per cell type per iteration (default: 100)",
            )
    parser.add_argument(
            "--n-iter", type=int, default=100,
            help="Stability iterations (default: 100)",
            )
    parser.add_argument(
            "--seed", type=int, default=20260904,
            help="Random seed, recorded in the summary for reproducibility",
            )
    parser.add_argument(
            "--aggregators", default=DEFAULT_AGGREGATORS,
            help=(
                f"Comma-separated marker aggregators (default: {DEFAULT_AGGREGATORS}). "
                "The first is primary and drives the stability loop; the rest are "
                "sensitivity checks."
                ),
            )
    parser.add_argument(
            "--expression", default="tpm", choices=("tpm", "tmm"),
            help=(
                "Which linear expression matrix to deconvolve (default: tpm). "
                "'tmm' uses edgeR's trimmed-mean scaling, which "
                "is robust to a minority of genes dominating the library -- worth "
                "comparing where a few transcripts hold much of the budget."
                ),
            )
    parser.add_argument(
            "--label", default=None, metavar="NAME",
            help=(
                "Suffix for this run's output files, e.g. --label 5type writes "
                "05_cell_proportions_5type.csv. Lets alternative cell-type sets "
                "coexist in one --output directory. Note params.json still records "
                "only the most recent run; each run's own 05_dsa_summary_<label>.json "
                "is the per-run provenance record."
                ),
            )
    parser.add_argument(
            "--stability-min", type=float, default=0.8,
            help=(
                "Median between-run correlation below which a cell type is reported as "
                "unstable (default: 0.8)"
                ),
            )
    parser.add_argument(
            "--proportions", default="auto",
            choices=("auto", "ensemble", "representative"),
            help=(
                "Which estimate becomes the reported proportions. 'ensemble' averages "
                "over all --n-iter marker draws and records the between-iteration SD; "
                "'representative' uses the single best-correlated draw. 'auto' (default) "
                "takes the ensemble when every retained cell type clears "
                "--stability-min and falls back to representative with a warning "
                "otherwise -- averaging over draws that disagree hides the disagreement "
                "rather than resolving it."
                ),
            )
    add_common_arguments(parser)


def run(args):
    output_dir = Path(args.output)
    for needed in ("02_expression.rds", "04_marker_coherence.csv"):
        if not (output_dir / needed).is_file():
            raise SystemExit(
                    f"tarva analysispipe: {output_dir / needed} not found. Run the "
                    f"expression and coherence steps into {output_dir} first."
                    )
    run_r_step(
            args,
            step="dsa",
            script_name="05_dsa.R",
            extra={
                "cell_types":    _csv_list(args.cell_types) if args.cell_types else None,
                "top_n":         args.top_n,
                "n_markers":     args.n_markers,
                "n_iter":        args.n_iter,
                "seed":          args.seed,
                "aggregators":   _csv_list(args.aggregators),
                "stability_min": args.stability_min,
                "proportions":   args.proportions,
                "expression":    args.expression,
                "label":         args.label,
                },
            )
