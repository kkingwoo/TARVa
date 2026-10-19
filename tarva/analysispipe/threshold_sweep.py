"""
analysispipe threshold-sweep: test how far the condition-enriched gene subsets
move when the enrichment threshold moves.

A gene joins a condition's subset when that condition holds at least 75% of
its alternate-allele calls; this repeats the assignment at 0.60, 0.70 and 0.80
and tests the four assignments for independence. Nothing is recomputed from
`call`: the counts do not depend on the threshold, only the label does, so it
reads gene_class's stored n_calls_case / n_calls_ref and re-derives the label
through gene_tests.assign_groups rather than with a second copy of it.

A NOTE ON THE TEST. The contingency table is thresholds x groups holding gene
counts, so every gene appears once per threshold: these are repeated measures
on one set of genes, not four independent samples, and chi-square's
independence assumption does not hold. Read the result as a descriptive
measure of how much the assignment moves, not as a hypothesis test. One test,
so BH adjustment is the identity and the raw p-value stands.

    `tarva analysispipe threshold-sweep --db <TARVa.db> [--output <dir>]`

Input:  a database with `gene-tests` done, which is where the per-gene counts
        and the threshold it used are read from.
Output: threshold_sweep.csv (Table S3) when --output is given, and
        threshold_sweep.* keys in run_meta.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import chi2_contingency

from tarva.analysispipe import gene_tests
from tarva.analysispipe import significance
from tarva.dbpipe import schema

# The three alternatives compared against the threshold in use. Not a flag:
# these three are fixed, and the fourth comes from whatever threshold
# gene-tests actually used.
SWEEP_THRESHOLDS = (0.60, 0.70, 0.80)

SWEEP_CSV = "threshold_sweep.csv"


def add_arguments(parser):
    parser.add_argument("--db", required=True,
                        help="Database with analysispipe gene-tests done")
    parser.add_argument("--output", default=None, metavar="DIR",
                        help="Directory for threshold_sweep.csv (Table S3). Created if "
                             "absent. Without it the result is reported to stdout and "
                             "stored in run_meta only")


def load_counts(con):
    """
    gene_class's per-gene alternate-allele counts, with the condition labels
    and the threshold gene-tests used -- all three from run_meta, so the sweep
    never has to re-derive which label was the case arm.
    """
    schema.require_table(con, "gene_class", "tarva analysispipe gene-tests")
    query = """
        SELECT ref_id, n_calls_case, n_calls_ref
        FROM gene_class
    """
    frame = pd.read_sql_query(query, con)
    case = schema.get_meta(con, "gene_tests.case_condition")
    reference = schema.get_meta(con, "gene_tests.reference_condition")
    baseline = schema.get_meta(con, "gene_tests.group_threshold")
    if not case or not reference or baseline is None:
        raise SystemExit(
                "tarva analysispipe threshold-sweep: run_meta holds no record of "
                "gene-tests' condition labels and group threshold, so the stored "
                "counts cannot be read. Re-run `tarva analysispipe gene-tests`.")
    return frame, str(case), str(reference), float(baseline)


def sweep(frame, case, reference, thresholds):
    """
    One row per (threshold, group): how many genes land there.

    Long form, because a cross-tab, a chi-square and a ggplot all read it
    without reshaping. Every group appears at every threshold, zero included:
    a missing cell and an empty cell are the same thing to chi-square, but
    only one of them can be plotted.
    """
    groups = [f"{case}_enriched", f"{reference}_enriched", "Shared"]
    records = []
    for threshold in thresholds:
        assigned = pd.Series(gene_tests.assign_groups(
                frame["n_calls_case"], frame["n_calls_ref"],
                case, reference, threshold))
        counts = assigned.value_counts()
        for group in groups:
            records.append({"threshold": threshold, "gene_group": group,
                            "n_genes": int(counts.get(group, 0))})
    return pd.DataFrame.from_records(records)


def test_independence(sweeps):
    """
    Chi-square of independence over the thresholds x groups table, with
    Cramer's V = sqrt(X2 / (n * (min(rows, cols) - 1))).

    n is every cell summed -- each gene once per threshold. A group empty at
    every threshold is dropped: it contributes nothing and chi-square cannot
    take a zero margin. Returns NaN statistics when that leaves fewer than two
    rows or columns, which is not an error but a cohort with no sensitivity to
    measure.
    """
    table = sweeps.pivot(index="threshold", columns="gene_group",
                         values="n_genes").fillna(0)
    table = table.loc[:, table.sum() > 0]
    table = table.loc[table.sum(axis=1) > 0, :]
    if min(table.shape) < 2:
        return table, np.nan, np.nan, 0, np.nan

    chi2, p_value, dof, _ = chi2_contingency(table.to_numpy())
    n = int(table.to_numpy().sum())
    denominator = n * (min(table.shape) - 1)
    cramers_v = float(np.sqrt(chi2 / denominator)) if denominator else np.nan
    return table, float(chi2), float(p_value), int(dof), cramers_v


def run(args):
    db_path = Path(args.db)
    if not db_path.is_file():
        raise SystemExit(
                f"tarva analysispipe threshold-sweep: database not found: {db_path}")

    con = schema.connect_for_build(db_path)
    try:
        frame, case, reference, baseline = load_counts(con)
        thresholds = sorted(set(SWEEP_THRESHOLDS) | {baseline})
        schema.report(f"thresholds: "
                      f"{', '.join(f'{t:.0%}' for t in thresholds)} "
                      f"(gene-tests assigned at {baseline:.0%})")
        schema.report(f"genes: {len(frame):,} from gene_class")

        # Carried forward rather than re-derived: if gene_class was built on
        # the all-genes fallback, so is this sweep, and that has to stay
        # attached to the number.
        filtered = schema.get_meta(con, "gene_tests.subset_from_global_set")
        if filtered is not None and not int(filtered):
            schema.warn("threshold sweep: gene_class was built over ALL genes because "
                        "step 1a's global set was empty, so this measures the "
                        "assignment's sensitivity on an unfiltered set. A mechanical "
                        "check of the code, not a filtered-set result.")

        sweeps = sweep(frame, case, reference, thresholds)
        table, chi2, p_value, dof, cramers_v = test_independence(sweeps)

        for threshold in table.index:
            assignments = table.loc[threshold]
            schema.report(f"  {threshold:.0%}  " + "   ".join(
                    f"{group} {int(n_genes):,}"
                    for group, n_genes in assignments.items()))

        if np.isnan(chi2):
            schema.warn("threshold sweep: the table has fewer than two non-empty rows "
                        "or columns, so chi-square is undefined. Every gene fell in "
                        "one group at every threshold -- the counts above are still "
                        "reportable, but there is no sensitivity to test.")
        else:
            schema.report(f"chi-square: X2 = {chi2:.4g}, df = {dof}, "
                          f"p = {significance.format_p(p_value)}, "
                          f"Cramer's V = {cramers_v:.3f}, over n = "
                          f"{int(table.to_numpy().sum()):,} gene-by-threshold "
                          f"observation(s)")

        schema.set_meta(con, "threshold_sweep", {
            "threshold_sweep.thresholds": ",".join(f"{t:g}" for t in thresholds),
            "threshold_sweep.baseline_threshold": baseline,
            "threshold_sweep.n_genes": len(frame),
            "threshold_sweep.chi2": chi2,
            "threshold_sweep.dof": dof,
            "threshold_sweep.p_value": p_value,
            "threshold_sweep.cramers_v": cramers_v,
            "threshold_sweep.output_dir": (str(Path(args.output).resolve())
                                           if args.output else ""),
            })
    finally:
        con.close()

    if args.output:
        directory = Path(args.output)
        directory.mkdir(parents=True, exist_ok=True)
        # The single test rides on every row. Table S3 is the counts table
        # with the statistic in its caption, and one file is what an .Rmd can
        # read without being told where the other half lives.
        out = sweeps.copy()
        out["chi2"] = chi2
        out["dof"] = dof
        out["p_value"] = p_value
        out["cramers_v"] = cramers_v
        gene_tests.write_csv(out, directory, SWEEP_CSV)

    schema.report("")
    schema.report("threshold-sweep complete")
    schema.report(f"  genes                     {len(frame):,}")
    schema.report(f"  thresholds                {len(thresholds)}")
    if not np.isnan(chi2):
        schema.report(f"  X2                        {chi2:.4g} (df {dof})")
        schema.report(f"  p                         {significance.format_p(p_value)}")
        schema.report(f"  Cramer's V                {cramers_v:.3f}")
    schema.report(f"Wrote {db_path}")
    return 0
