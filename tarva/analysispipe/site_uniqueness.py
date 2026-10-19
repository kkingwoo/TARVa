"""
analysispipe site-uniqueness: are uniquely modified positions distributed
differently between the two conditions, per variant call-type and per gene
subset?

    unique   a position called in EXACTLY ONE individual, counted across both
             conditions together
    common   a position called in more than one individual, whether within or
             across conditions

Uniqueness is a property of the POSITION, decided once over the whole cohort,
not per condition. A common position called in both conditions is therefore
counted in both conditions' `common` cell, while a unique position belongs to
exactly one condition's `unique` cell by construction.

Then per (gene subset, call type), a 2x2 of condition against uniqueness,
counting DISTINCT POSITIONS:

                  unique            common
    case     n_case_unique     n_case_common
    ref      n_ref_unique      n_ref_common

Two-sided Fisher's exact on those counts, then Benjamini-Hochberg across every
(subset, call type) the step computed.

THE GENE SUBSET IS edit_dist, NOT gene_group -- the per-gene majority-rule
distribution label rather than the 75% subset assignment. Both name two of
their three categories `*_enriched`, so they are easy to confuse.

ODDS RATIO AND INTERVAL. The ratio is
(case_unique * ref_common) / (case_common * ref_unique), so above 1 means
unique sites are relatively more likely in the case condition. It is computed
on counts with Haldane-Anscombe applied -- a zero cell is replaced by 0.5 --
because a zero makes the ratio infinite or undefined. The interval is Woolf's:
exp(log(OR) +/- 1.96 * sqrt(sum of 1/cell)). The p-value, by contrast, is
computed on the RAW counts: Fisher's exact needs no continuity or zero
correction, and correcting its input would change a p-value that was already
valid.

    `tarva analysispipe site-uniqueness --db <TARVa.db>
        [--min-depth 20] [--output <dir>] [--overwrite]`

Input:  a database with `gene-tests` done, for the gene subsets and the
        condition labels.
Output: the uniqueness_test table, site_uniqueness.csv (Figure 5, Table S1)
        when --output is given, and site_uniqueness.* keys in run_meta.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import fisher_exact

from tarva.analysispipe import gene_tests
from tarva.analysispipe import significance
from tarva.dbpipe import schema

UNIQUENESS_TABLES = ("uniqueness_test",)

UNIQUENESS_CSV = "site_uniqueness.csv"

UNIQUE, COMMON = "unique", "common"

# Haldane-Anscombe: a zero cell becomes this, for the ratio and interval only.
HALDANE = 0.5

# The normal quantile for a 95% interval. 1.96, not 1.959964.
Z_95 = 1.96

DEFAULT_MIN_DEPTH = 20
ALPHA = 0.05


def add_arguments(parser):
    parser.add_argument("--db", required=True,
                        help="Database with gene-tests done")
    parser.add_argument("--min-depth", type=int, default=DEFAULT_MIN_DEPTH,
                        help=f"A call counts toward uniqueness only at or above this "
                             f"depth (default: {DEFAULT_MIN_DEPTH}). Should match the "
                             f"--min-depth gene-tests used")
    parser.add_argument("--output", default=None, metavar="DIR",
                        help="Directory for site_uniqueness.csv (Figure 5, Table S1). "
                             "Created if absent. Without it the result goes to the "
                             "database and stdout only")
    parser.add_argument("--overwrite", action="store_true",
                        help="Drop and recompute uniqueness_test instead of failing")


#----------------------------------------------------
# inputs
#----------------------------------------------------

def load_sites(con, min_depth):
    """
    One row per (position, individual): the gene's edit_dist, the call type,
    and the condition.

    Joined to gene_class rather than gene, which restricts the result to the
    gene set steps 2a and 2b used and supplies edit_dist in the same pass.
    """
    schema.require_table(con, "gene_class", "tarva analysispipe gene-tests")
    schema.require_table(con, "call", "tarva dbpipe load-calls")
    table = schema.INDIVIDUAL_TABLE
    id_column, condition_column = schema.individual_columns(con)
    query = f"""
        SELECT s.chrom          AS chrom,
               s.pos            AS pos,
               g.edit_dist      AS edit_dist,
               c.a1type         AS a1type,
               c.individualID   AS individualID,
               p."{condition_column}" AS condition
        FROM call AS c
        JOIN site AS s ON s.chrom = c.chrom AND s.pos = c.pos
        JOIN gene_class AS g ON g.ref_id = s.ref_id
        JOIN {table} AS p ON p."{id_column}" = c.individualID
        WHERE c.a1type IS NOT NULL
          AND c.dp >= ?
    """
    frame = pd.read_sql_query(query, con, params=(int(min_depth),))
    if frame.empty:
        raise SystemExit(
                f"tarva analysispipe site-uniqueness: no call in the gene_class gene "
                f"set has a call type and dp >= {min_depth}. Nothing to classify.")
    return frame


def classify_uniqueness(frame):
    """
    Label each row by whether its POSITION is carried by exactly one
    individual across the whole cohort.

    Counted over both conditions together, which is what makes `unique` mean
    "nobody else has it" rather than "nobody else in this condition has it".
    A position's label is a property of the position, so every row at that
    position gets the same one.
    """
    carriers = frame.groupby(["chrom", "pos"])["individualID"].transform("nunique")
    frame = frame.copy()
    frame["n_carriers"] = carriers
    frame["uniqueness"] = np.where(carriers == 1, UNIQUE, COMMON)
    return frame


def count_positions(frame):
    """
    DISTINCT positions per (edit_dist, a1type, condition, uniqueness).

    De-duplicated before the count, so a position called by several
    individuals of the same condition counts once for that condition.
    """
    keys = ["edit_dist", "a1type", "condition", "uniqueness"]
    distinct = frame.drop_duplicates(keys + ["chrom", "pos"])
    return distinct.groupby(keys).size().rename("n_positions").reset_index()


#----------------------------------------------------
# the test
#----------------------------------------------------

def odds_ratio(case_unique, case_common, ref_unique, ref_common):
    """
    (OR, ci_low, ci_high) on Haldane-Anscombe-corrected counts.

    Returns NaN throughout when a whole row or column of the table is zero
    even after correction is impossible to rescue -- there is no ratio to
    report when one condition contributed no positions at all.
    """
    cells = [float(n) for n in (case_unique, case_common, ref_unique, ref_common)]
    if not any(cells):
        return np.nan, np.nan, np.nan
    corrected = [HALDANE if cell == 0 else cell for cell in cells]
    a, b, c, d = corrected
    ratio = (a * d) / (b * c)
    se = np.sqrt(sum(1.0 / cell for cell in corrected))
    return (ratio,
            float(np.exp(np.log(ratio) - Z_95 * se)),
            float(np.exp(np.log(ratio) + Z_95 * se)))


def run_tests(counts, case, reference):
    """
    One row per (edit_dist, a1type): the 2x2, Fisher's exact, OR and interval.

    Every (subset, call type) present in the data gets a row, with absent
    cells as 0 rather than missing, so a call type seen only in one condition
    still appears -- it is a result, not a gap.
    """
    wide = counts.pivot_table(
            index=["edit_dist", "a1type"],
            columns=["condition", "uniqueness"],
            values="n_positions",
            fill_value=0,
            aggfunc="sum")

    records = []
    for (edit_dist, a1type), row in wide.iterrows():
        cells = {}
        for label, prefix in ((case, "case"), (reference, "ref")):
            for kind in (UNIQUE, COMMON):
                cells[f"n_{prefix}_{kind}"] = int(row.get((label, kind), 0))

        case_unique = cells["n_case_unique"]
        case_common = cells["n_case_common"]
        ref_unique = cells["n_ref_unique"]
        ref_common = cells["n_ref_common"]

        # Fisher's exact on the RAW counts; it handles zeros exactly.
        _, p_value = fisher_exact(
                [[case_unique, case_common], [ref_unique, ref_common]],
                alternative="two-sided")
        ratio, ci_low, ci_high = odds_ratio(
                case_unique, case_common, ref_unique, ref_common)

        records.append(dict(
            edit_dist=edit_dist, a1type=a1type,
            odds_ratio=ratio, ci_low=ci_low, ci_high=ci_high,
            p_value=p_value, **cells))

    tests = pd.DataFrame.from_records(records)
    if tests.empty:
        return tests
    tests["adj_pvalue"] = significance.benjamini_hochberg(
            tests["p_value"].to_numpy())
    return tests


#----------------------------------------------------
# writing
#----------------------------------------------------

def check_reloadable(con, overwrite):
    counts = {table: con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
              for table in UNIQUENESS_TABLES}
    if any(counts.values()) and not overwrite:
        raise SystemExit(
                f"tarva analysispipe site-uniqueness: "
                f"{', '.join(f'{t} holds {n:,} row(s)' for t, n in counts.items())}. "
                f"Pass --overwrite to recompute.")


def prepare_tables(con):
    counts = {table: con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
              for table in UNIQUENESS_TABLES}
    if any(counts.values()):
        for table in UNIQUENESS_TABLES:
            con.execute(f"DROP TABLE IF EXISTS {table}")
        con.commit()
        schema.create_schema(con, UNIQUENESS_TABLES)
        schema.report(f"--overwrite: dropped and recreated "
                      f"{', '.join(UNIQUENESS_TABLES)}")


#----------------------------------------------------
# entry point
#----------------------------------------------------

def run(args):
    db_path = Path(args.db)
    if not db_path.is_file():
        raise SystemExit(
                f"tarva analysispipe site-uniqueness: database not found: {db_path}")

    con = schema.connect_for_build(db_path)
    try:
        schema.create_schema(con)
        check_reloadable(con, args.overwrite)

        case, reference, sizes = gene_tests.load_conditions(con, None)
        schema.report(f"conditions: case={case} (n={sizes[case]}), "
                      f"reference={reference} (n={sizes[reference]})")

        filtered = schema.get_meta(con, "gene_tests.subset_from_global_set")
        if filtered is not None and not int(filtered):
            schema.warn("site-uniqueness: gene_class was built over ALL genes because "
                        "step 1a's global set was empty, so these counts describe the "
                        "unfiltered gene set. run_meta records which set was used.")

        frame = classify_uniqueness(load_sites(con, args.min_depth))
        positions = frame.drop_duplicates(["chrom", "pos"])
        n_unique = int((positions["uniqueness"] == UNIQUE).sum())
        n_common = int((positions["uniqueness"] == COMMON).sum())
        schema.report(f"positions: {len(positions):,} distinct, {n_unique:,} unique "
                      f"(one carrier), {n_common:,} common")
        schema.report(f"carriers per position: median "
                      f"{positions['n_carriers'].median():.0f}, max "
                      f"{int(positions['n_carriers'].max())}")

        counts = count_positions(frame)
        tests = run_tests(counts, case, reference)
        schema.report(f"tests: {len(tests):,} (subset, call type) over "
                      f"{tests['edit_dist'].nunique()} subset(s) and "
                      f"{tests['a1type'].nunique()} call type(s)")

        basis, n_significant = significance.choose_basis(
                tests, ALPHA, significance.REPORTING, "site uniqueness",
                downstream="Figure 5's bolded rows")
        # Everything below reports on the basis the cascade chose, not on
        # adj_pvalue regardless. Reporting a nominal count under an "adjusted"
        # label is how a fallback becomes invisible.
        basis_label = {significance.ADJUSTED: "adjusted",
                       significance.RAW: "nominal"}.get(basis, "no qualifying basis")
        tests["basis"] = basis or ""

        for edit_dist, block in tests.groupby("edit_dist", sort=True):
            significant = (block.iloc[0:0] if basis is None
                           else block.loc[block[basis] < ALPHA])
            schema.report(f"  {edit_dist}: {len(block):,} call type(s), "
                          f"{len(significant):,} at {basis_label} p < {ALPHA}")
            for row in significant.sort_values("p_value").itertuples():
                direction = (case if row.odds_ratio > 1 else reference)
                schema.report(f"    {row.a1type:<8} OR {row.odds_ratio:>8.3f} "
                              f"[{row.ci_low:.3f}, {row.ci_high:.3f}]  "
                              f"unique more likely in {direction}")

        prepare_tables(con)
        tests.reindex(columns=[
            "edit_dist", "a1type", "n_case_unique", "n_case_common",
            "n_ref_unique", "n_ref_common", "odds_ratio", "ci_low", "ci_high",
            "p_value", "adj_pvalue", "basis"]).to_sql(
                "uniqueness_test", con, if_exists="append", index=False,
                chunksize=gene_tests.WRITE_CHUNK)
        con.commit()
        schema.report(f"uniqueness_test: {len(tests):,} row(s) written")

        schema.create_indexes(con, UNIQUENESS_TABLES)
        schema.set_meta(con, "site_uniqueness", {
            "site_uniqueness.case_condition": case,
            "site_uniqueness.reference_condition": reference,
            "site_uniqueness.min_depth": args.min_depth,
            "site_uniqueness.alpha": ALPHA,
            "site_uniqueness.haldane_anscombe": HALDANE,
            "site_uniqueness.n_positions": len(positions),
            "site_uniqueness.n_unique": n_unique,
            "site_uniqueness.n_common": n_common,
            "site_uniqueness.n_tests": len(tests),
            "site_uniqueness.basis": basis or "none",
            "site_uniqueness.n_significant": n_significant,
            "site_uniqueness.output_dir": (str(Path(args.output).resolve())
                                           if args.output else ""),
            })
    finally:
        con.close()

    if args.output:
        directory = Path(args.output)
        directory.mkdir(parents=True, exist_ok=True)
        # The condition labels ride along so the forest plot can name its
        # sections without being told which group was the case arm.
        out = tests.copy()
        out["case_condition"] = case
        out["ref_condition"] = reference
        gene_tests.write_csv(out, directory, UNIQUENESS_CSV)

    schema.report("")
    schema.report("site-uniqueness complete")
    schema.report(f"  distinct positions        {len(positions):,}")
    schema.report(f"  unique                    {n_unique:,}")
    schema.report(f"  common                    {n_common:,}")
    schema.report(f"  (subset, call type) tests {len(tests):,}")
    schema.report(f"  significant               {n_significant:,} "
                  f"({basis_label} p < {ALPHA})")
    schema.report(f"Wrote {db_path}")
    return 0
