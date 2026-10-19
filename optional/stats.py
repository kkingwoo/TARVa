"""
analysispipe stats: group comparison, neuropathology correlation and covariate
    adjustment for the estimated cell-type proportions.

    `tarva analysispipe stats --db <tarva.db> --output <results_dir>
        [--label 5type_tmm] [--case-group AD]
        [--continuous braaksc] [--covariates age_death,msex,pmi]
        [--covariate-file <pmi.csv>] [--covariate-file-key individualID]
        [--aggregator median]
        [--hpc] [--r-module R/4.5.1] [--r-libs PATH] [--directives-dict "{...}"]`

    Input:  the proportions `dsa` wrote into the same --output under the same
            --label, and the database for the cohort's covariates.
    Output: 06_covariates.csv, 06_stats_binary.csv, 06_stats_continuous.csv,
            06_stats_adjusted.csv, 06_stats_summary.json,
            06_proportions_by_group.png

    Three analyses, because the binary contrast alone understates what the data
    can say:

      binary     -- medians and IQRs per group, case:control ratio, Mann-Whitney
                    with BH correction across cell types, Hodges-Lehmann shift
                    with CI, and Cliff's delta.

      continuous -- Spearman of each proportion against neuropathology. Better
                    powered than the dichotomy whenever the groups overlap on
                    pathology, and it answers the more direct question: does
                    composition track the disease process at all?

      adjusted   -- beta regression, since proportions are bounded in (0,1).
                    Reported beside the unadjusted result, not instead of it.

    A covariate perfectly collinear with the grouping is dropped from the
    adjusted model with a warning -- it carries no within-group variation and
    would leave the model unidentifiable.

    Effect sizes and confidence intervals are the output that matters. At
    modest n a non-significant result must read as "consistent with a small
    change in either direction", never as "no change".
"""

from pathlib import Path

from tarva.analysispipe import db as tarva_db
from tarva.analysispipe.rstep import add_common_arguments, run_r_step

# braaksc is deliberately absent. If the causal path is diagnosis -> tangle
# pathology -> cell loss, conditioning on Braak removes the effect being
# measured. It is fitted separately as --mediator-covariates and labelled.
DEFAULT_COVARIATES = "age_death,msex,pmi"
DEFAULT_MEDIATORS = "braaksc"


def _csv_list(raw):
    return [v.strip() for v in str(raw).split(",") if v.strip()]


def add_arguments(parser):
    parser.add_argument(
            "--db", required=True,
            help="TARVa SQLite database -- supplies the covariates",
            )
    parser.add_argument(
            "--table", default=tarva_db.DEFAULT_TABLE,
            help=f"Database table holding sample info (default: {tarva_db.DEFAULT_TABLE})",
            )
    parser.add_argument(
            "--label", default=None, metavar="NAME",
            help="Label of the dsa run to analyse, matching that run's --label",
            )
    parser.add_argument(
            "--case-group", default=None,
            help=(
                "Value of `con` treated as the case group. Defaults to the first "
                "sorted value, which is printed -- check it, since it sets the "
                "direction of every ratio and shift."
                ),
            )
    parser.add_argument(
            "--continuous", default="braaksc",
            help=(
                "Comma-separated continuous variables to correlate proportions against "
                "(default: braaksc). A variable that perfectly separates the groups "
                "carries no extra information and is better used as the grouping itself."
                ),
            )
    parser.add_argument(
            "--covariates", default=DEFAULT_COVARIATES,
            help=(
                f"Comma-separated covariates for the PRIMARY adjusted model "
                f"(default: {DEFAULT_COVARIATES}). Should contain no variable on the "
                f"causal path from diagnosis to composition."
                ),
            )
    parser.add_argument(
            "--mediator-covariates", default=DEFAULT_MEDIATORS,
            help=(
                f"Comma-separated variables fitted as a SEPARATE, labelled sensitivity "
                f"model (default: {DEFAULT_MEDIATORS}). These are plausibly on the causal "
                f"path -- conditioning on them removes part of the effect being measured, "
                f"so the result is reported with that caveat rather than as the adjusted "
                f"estimate. Pass an empty string to skip."
                ),
            )
    parser.add_argument(
            "--aggregator", default="median",
            help=(
                "Aggregator treated as primary in the plot and summary (default: median, "
                "per Wang et al.). All aggregators present are analysed regardless."
                ),
            )
    tarva_db.add_covariate_file_arguments(parser)
    add_common_arguments(parser)


def run(args):
    output_dir = Path(args.output)
    sfx = f"_{args.label}" if args.label else ""
    props = output_dir / f"05_cell_proportions{sfx}.csv"
    if not props.is_file():
        raise SystemExit(
                f"tarva analysispipe: {props} not found. Run `tarva analysispipe dsa "
                f"--output {output_dir}"
                + (f" --label {args.label}" if args.label else "")
                + "` first."
                )

    covariates = tarva_db.load_covariates(args.db, args.table)
    if args.covariate_file:
        covariates = tarva_db.merge_covariate_file(
                covariates, args.covariate_file,
                key=args.covariate_file_key,
                columns=_csv_list(args.covariate_file_columns) if args.covariate_file_columns else None,
                )
    covariates.to_csv(output_dir / "06_covariates.csv", index=False)
    print(f"Loaded covariates for {len(covariates)} samples: "
          + ", ".join(c for c in covariates.columns if c not in ("rids", "con")))

    run_r_step(
            args,
            step="stats",
            script_name="06_stats.R",
            extra={
                "label":       args.label,
                "case_group":  args.case_group,
                "continuous":  _csv_list(args.continuous),
                "covariates":  _csv_list(args.covariates),
                "mediators":   _csv_list(args.mediator_covariates) if args.mediator_covariates else [],
                "aggregator":  args.aggregator,
                },
            job_name=f"tarva_stats{sfx}",
            )
