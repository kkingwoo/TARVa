"""
analysispipe confounding: does estimated cell composition account for the
    editing findings?

    `tarva analysispipe confounding --db <tarva.db> --gene-groups <TableS4.csv>
        --output <results_dir> [--label 5type_tmm] [--case-group AD]
        [--aggregator median] [--drop-cell-type neu] [--seed 20260904]
        [--site-genes GENE1,GENE2,...]
        [--hpc] [--r-module R/4.5.1] [--r-libs PATH] [--directives-dict "{...}"]`

    Input:  the database, a CSV of gene class labels, the proportions `dsa`
            wrote into the same --output under the same --label, and
            08_gene_group_resolved.csv from `check_gene_groups`.
    Output: 08_burden.csv, 08_gene_sites.csv, 08_site_level.csv,
            08_burden_vs_composition.csv, 08_gene_models.csv,
            08_class_equivalence.csv, 08_class_comparison.csv,
            08_call_models.csv, 08_call_sensitivity.csv, 08_depth_models.csv,
            08_site_models.csv, 08_confounding_summary.json,
            08_confounding.png

    The class labels in --gene-groups are taken AS GIVEN. Nothing here
    re-derives them, and no enrichment threshold appears anywhere in this
    step: a label re-invented from a single fraction-of-sites cut will not
    match the one the gene-level analysis assigned. The equivalent question --
    does a gene's apparent effect survive accounting for how many calls an
    individual yields overall -- is asked as a model covariate instead.

    Estimating composition answers only half the question: confounding
    requires the confounder to differ between groups AND to associate with the
    outcome. The components below address the mechanism at the level where it
    would operate rather than at its traces.

    THE MECHANISM. Calls are made only where dp >= --dp-threshold. So:
    composition shifts -> expression of a cell type's genes shifts -> read depth
    at their positions shifts -> fewer or more positions clear the threshold ->
    the gene appears condition-enriched. Presence/absence findings are the most
    exposed, because a detection-sensitivity shift can manufacture them outright.

    WHAT THIS DATA CAN SETTLE. dp is recorded only where a call was made, so a
    (gene, sample) pair with no call is either a real negative or a position
    that was never testable. Nothing here can tell those apart -- that needs the
    VCFs or BAMs. Rather than impute which zeros were untestable, the models
    adjust for the OBSERVED determinants of depth: gene expression, which exists
    for every pair, and library size. Every model is fitted within a gene, so
    sample-constant coverage bias is absorbed; what escapes is sample-by-gene
    interaction in coverage bias. There is no basis here for a single
    attributable-fraction figure and none is produced.

    COMPONENTS.
      8a  Per-individual call burden against condition and against each
          estimated proportion, marginally and partialled on diagnosis. The one
          component with no missing-data problem: burden is a property of the
          individual and is not conditional on any gene being testable. If
          burden is unrelated to composition the mechanism has no foothold.

      8b  Per-gene expression regressed on composition, adjusted R2 and one
          coefficient per cell type, compared across the enriched/shared
          classes with Shared as the built-in null. A proxy for confounding of
          the calls, not a measurement of it -- 8f is that.

      8c  Per-cell-type coefficients compared across classes, one-sided in the
          direction the artefact hypothesis predicts.

      8f  THE DIRECT TEST. Per gene, whether an individual carries any call,
          modelled on diagnosis alone, then plus composition, then plus
          expression and library size. The movement in the diagnosis
          coefficient is the magnitude of confounding for that gene. Reported
          with a sensitivity table across minimum-expression cut-offs so no
          single risk-set definition is load-bearing. Firth's penalised
          likelihood is used when logistf is installed, since genes called in
          one group only separate perfectly.

      8g  Observed read depth against composition, among called positions only.
          A lower bound on the composition-to-depth arrow; feeds nothing.

      8d  Site-level editing proportions, quasibinomial, naive versus
          composition-adjusted, across every labelled gene.

      8e  Named genes located within 8b's and 8f's distributions.

    Proportions sum to 1, so one cell type is dropped from every model --
    including all of them is perfectly collinear.
"""

from pathlib import Path

from tarva.analysispipe import db as tarva_db
from tarva.analysispipe.rstep import add_common_arguments, run_r_step


def _csv_list(raw):
    return [v.strip() for v in str(raw).split(",") if v.strip()]


def add_arguments(parser):
    parser.add_argument("--db", required=True, help="TARVa SQLite database")
    parser.add_argument(
            "--table", default=tarva_db.DEFAULT_TABLE,
            help=f"Database table holding the editing calls (default: {tarva_db.DEFAULT_TABLE})",
            )
    parser.add_argument(
            "--gene-groups", required=True,
            help="CSV of globally-labelled genes, used to locate 08_gene_group_resolved.csv",
            )
    parser.add_argument(
            "--label", default=None, metavar="NAME",
            help="Label of the dsa run to use, matching that run's --label",
            )
    parser.add_argument(
            "--case-group", default=None,
            help="Value of `con` treated as the case group (default: first sorted, printed)",
            )
    parser.add_argument(
            "--aggregator", default="median",
            help="Proportion aggregator to use (default: median)",
            )
    parser.add_argument(
            "--expression", default=None, choices=("tpm", "tmm"),
            help=(
                "Expression matrix to model against. Defaults to whatever the labelled "
                "dsa run used, read from its own summary -- modelling composition "
                "against a differently-normalised matrix than the proportions came from "
                "puts a per-sample scale difference into every coefficient. The "
                "unfiltered variant is used either way, since it carries the labelled "
                "genes filterByExpr dropped."
                ),
            )
    parser.add_argument(
            "--drop-cell-type", default=None,
            help=(
                "Cell type omitted from every model to avoid perfect collinearity "
                "(proportions sum to 1). Default: the most abundant, which is the "
                "conventional choice."
                ),
            )
    parser.add_argument(
            "--seed", type=int, default=20260904,
            help="Random seed for the equivalence bootstrap, recorded in the summary",
            )
    parser.add_argument(
            "--min-n", type=int, default=20,
            help=(
                "Minimum samples for a site-level model to be treated as estimable "
                "(default: 20). The composition-adjusted model fits six parameters, so "
                "below this the residual degrees of freedom leave the coefficient "
                "uninterpretable. Positions under the threshold are still modelled and "
                "written to the CSV, but are excluded from the BH family and reported "
                "separately rather than counted as attenuation."
                ),
            )
    parser.add_argument(
            "--r2-equivalence-margin", type=float, default=0.10,
            help=(
                "Margin for the equivalence test on between-class differences in median "
                "composition R2 (default: 0.10). A bootstrap CI falling entirely within "
                "+/- this margin supports equivalence; a non-significant difference on "
                "its own does not. Pre-specify this and report it."
                ),
            )
    parser.add_argument(
            "--n-bootstrap", type=int, default=2000,
            help="Bootstrap resamples for the class-difference CI (default: 2000)",
            )
    parser.add_argument(
            "--site-genes", default=None,
            help=(
                "Comma-separated gene symbols printed in full in the per-gene detail "
                "report (8e). The adjusted models in 8d and 8f run on EVERY labelled "
                "gene regardless -- this only selects what is shown, not what is fitted."
                ),
            )
    parser.add_argument(
            "--min-carriers", type=int, default=5,
            help=(
                "Minimum individuals with, and without, a call in a gene for its "
                "gene-level presence model to be fitted (default: 5). Below this the "
                "logistic coefficient is driven by a handful of people and measures "
                "nothing."
                ),
            )
    parser.add_argument(
            "--dp-threshold", type=int, default=20,
            help=(
                "Read depth at or above which a position was callable (default: 20). "
                "Recorded in the summary and used in the reported framing of the "
                "detection mechanism; set it to whatever the variant caller actually "
                "used."
                ),
            )
    parser.add_argument(
            "--case-enriched-group", default=None, metavar="LABEL",
            help=(
                "Value of the gene_group column marking case-enriched genes. Needed for "
                "the directional test in 8c, which is what the artefact hypothesis "
                "predicts; without it only the omnibus test is reported. Inferred from "
                "--case-group by name match when omitted."
                ),
            )
    parser.add_argument(
            "--ref-enriched-group", default=None, metavar="LABEL",
            help="Value of gene_group marking reference-enriched genes; see --case-enriched-group",
            )
    add_common_arguments(parser)


def run(args):
    output_dir = Path(args.output)
    sfx = f"_{args.label}" if args.label else ""

    resolved = output_dir / "08_gene_group_resolved.csv"
    if not resolved.is_file():
        raise SystemExit(
                f"tarva analysispipe: {resolved} not found. Run "
                f"`python -m tarva.analysispipe.check_gene_groups --gene-groups "
                f"{args.gene_groups} --output {output_dir}` first."
                )
    props = output_dir / f"05_cell_proportions{sfx}.csv"
    if not props.is_file():
        raise SystemExit(f"tarva analysispipe: {props} not found -- run the dsa step first.")

    import pandas as pd
    groups = pd.read_csv(resolved)
    genes = sorted(set(groups["gene"].astype(str)) | set(groups["ensgid"].astype(str)))

    burden = tarva_db.load_burden(args.db, args.table)
    burden.to_csv(output_dir / "08_burden.csv", index=False)
    print(f"Call burden for {len(burden)} samples: "
          f"median {int(burden['n_calls'].median())} calls, "
          f"{int(burden['n_positions'].median())} positions")

    sites = tarva_db.load_gene_sites(args.db, genes, args.table)
    sites.to_csv(output_dir / "08_gene_sites.csv", index=False)
    print(f"Loaded {len(sites)} (gene, sample, position) triples "
          f"across {sites['ref_base'].nunique()} genes")

    # Site-level rows for EVERY labelled gene, not only the named ones. The
    # composition-adjusted models in 8d and 8f are the only place a CN/AD effect
    # is re-estimated with composition in it; gating that on a hand-typed gene
    # list left the central question answered for a handful of genes and
    # unanswered for the rest.
    level = tarva_db.load_site_level(args.db, genes, args.table)
    level.to_csv(output_dir / "08_site_level.csv", index=False)
    print(f"Loaded {len(level)} site-level rows across "
          f"{level['ref_base'].nunique()} labelled gene(s)")

    site_genes = _csv_list(args.site_genes) if args.site_genes else []

    run_r_step(
            args,
            step="confounding",
            script_name="08_confounding.R",
            extra={
                "label":            args.label,
                "case_group":       args.case_group,
                "aggregator":       args.aggregator,
                # Nulled rather than omitted: params.json accumulates across steps,
                # so an --expression left by the last dsa run would otherwise take
                # precedence over the labelled run's own summary and silently select
                # the wrong matrix whenever more than one label exists.
                "expression":       args.expression,
                "drop_cell_type":   args.drop_cell_type,
                "seed":             args.seed,
                "min_n":            args.min_n,
                "min_carriers":     args.min_carriers,
                "dp_threshold":     args.dp_threshold,
                "r2_margin":        args.r2_equivalence_margin,
                "n_bootstrap":      args.n_bootstrap,
                "site_genes":       site_genes,
                "case_enriched_group": args.case_enriched_group,
                "ref_enriched_group":  args.ref_enriched_group,
                },
            job_name=f"tarva_confounding{sfx}",
            )
