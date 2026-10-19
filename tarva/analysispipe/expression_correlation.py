"""
analysispipe expression-correlation: for each gene, across the individuals
within a condition, is expression associated with variant load -- and does
that association differ between the two conditions?

    per gene, per condition:  Spearman(TPM, variant count) over individuals
    per gene:                 Fisher z on the difference of the two rho

Both quantities already exist in the database, so nothing is re-parsed:

  TPM     SUM(expression.TPM) over the gene's transcripts. `load-expression`
          already normalised FPKM over the whole quantified transcriptome
          before any filtering: all FPKM summed per individual, each gene's
          FPKM divided by that total, times 1e6.
  varcount COUNT of that individual's eligible calls in the gene, under the
          same eligibility every other step uses: resolved to a gene, a
          call type assigned, dp >= --min-depth.

EVERY GENE GETS A ROW FOR EVERY INDIVIDUAL. An individual with no calls in a
gene has a variant count of 0, and one with no quantified transcript has a TPM
of 0; both are measurements rather than gaps, so the correlation runs over the
whole condition.

The per-condition entry rule is that the condition must have at least
MIN_INDIVIDUALS individuals AND more than one distinct TPM AND more than one
distinct variant count. Spearman is undefined or unstable otherwise, so that
condition's rho and p are NULL -- a different statement from a correlation of
zero, and one that keeps the gene out of the FDR step rather than counting it
as non-significant.

NOTE FOR SMALL COHORTS. Because every individual is present, the arm size IS
the condition's cohort size. A condition with fewer than MIN_INDIVIDUALS
individuals therefore fails the rule for EVERY gene, and the step produces
nothing at all. That is a property of the cohort, not of the data.

    `tarva analysispipe expression-correlation --db <TARVa.db>
        [--min-depth 20] [--output <dir>] [--overwrite]`

Input:  a database with `gene-tests` done. The genes are those `gene_class`
        holds, so this describes the same gene set as the other steps, and
        `gene_group` comes along for the figure's facets.
Output: the expression_test table, tpm_variant_relationship.csv (Figure 4,
        Table S4) when --output is given, and expression_correlation.* keys in
        run_meta.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import norm, spearmanr

from tarva.analysispipe import gene_tests
from tarva.analysispipe import significance
from tarva.dbpipe import schema

EXPRESSION_TABLES = ("expression_test",)

RELATIONSHIP_CSV = "tpm_variant_relationship.csv"

# Not a flag: one value, applied everywhere.
#
# READ THIS BEFORE CHANGING IT. Because build_matrix zero-fills the full
# gene x individual cross product, len(arm) is ALWAYS the condition's cohort
# size, identical for every gene. So this is not a per-gene filter, however
# much the comparison in `correlate` looks like one: it is a gate on whether
# the step runs AT ALL. Raising it to n means the step refuses on any cohort
# with fewer than n individuals in a condition, for every gene at once.
MIN_INDIVIDUALS = 4

DEFAULT_MIN_DEPTH = 20

# The difference of correlations is the result here, not a filter, so it is
# reported at this alpha rather than gating anything.
DIFF_ALPHA = 0.05


def add_arguments(parser):
    parser.add_argument("--db", required=True,
                        help="Database with gene-tests done")
    parser.add_argument("--min-depth", type=int, default=DEFAULT_MIN_DEPTH,
                        help=f"A call counts toward a gene's variant load only at or "
                             f"above this depth (default: {DEFAULT_MIN_DEPTH}). Should "
                             f"match the --min-depth gene-tests used")
    parser.add_argument("--output", default=None, metavar="DIR",
                        help="Directory for tpm_variant_relationship.csv (Figure 4, "
                             "Table S4). Created if absent. Without it the result goes "
                             "to the database and stdout only")
    parser.add_argument("--overwrite", action="store_true",
                        help="Drop and recompute expression_test instead of failing")


#----------------------------------------------------
# inputs
#----------------------------------------------------

def load_gene_set(con):
    """
    The genes steps 2a and 2b used, with their subset labels.

    Read from gene_class rather than re-derived from gene_test, since that
    table already holds exactly the resolved gene set -- so this step cannot
    end up describing a different set of genes from the others.
    """
    schema.require_table(con, "gene_class", "tarva analysispipe gene-tests")
    query = """
        SELECT c.ref_id, c.gene_group, g.novel
        FROM gene_class AS c
        JOIN gene AS g ON g.ref_id = c.ref_id
    """
    return pd.read_sql_query(query, con)


def load_tpm(con, genes):
    """
    Per-gene TPM per individual: SUM over the gene's transcripts.

    Summed in SQL rather than in pandas because `expression` carries one row
    per transcript per individual and only the per-gene total is wanted.
    """
    schema.require_table(con, "expression", "tarva dbpipe load-expression")
    query = """
        SELECT e.individualID AS individualID,
               t.ref_id       AS ref_id,
               SUM(e.TPM)     AS tpm
        FROM expression AS e
        JOIN transcript AS t ON t.tid = e.tid
        WHERE t.ref_id IS NOT NULL
        GROUP BY e.individualID, t.ref_id
    """
    frame = pd.read_sql_query(query, con)
    return frame.loc[frame["ref_id"].isin(genes)]


def load_varcount(con, genes, min_depth):
    """
    Per-gene variant count per individual, under the pipeline's one
    eligibility rule.

    The same three conditions every other step applies -- resolved to a gene,
    a call type assigned, dp at or above --min-depth -- so "variant load"
    means the same thing here as everywhere else.
    """
    schema.require_table(con, "call", "tarva dbpipe load-calls")
    query = """
        SELECT c.individualID AS individualID,
               s.ref_id       AS ref_id,
               COUNT(*)       AS varcount
        FROM call AS c
        JOIN site AS s ON s.chrom = c.chrom AND s.pos = c.pos
        WHERE s.ref_id IS NOT NULL
          AND c.a1type IS NOT NULL
          AND c.dp >= ?
        GROUP BY c.individualID, s.ref_id
    """
    frame = pd.read_sql_query(query, con, params=(int(min_depth),))
    return frame.loc[frame["ref_id"].isin(genes)]


def build_matrix(genes, individuals, tpm, varcount, condition_of):
    """
    One row per (gene, individual): tpm, varcount, condition. Zero-filled.

    The full cross product is built first and the two measurements merged
    into it, so an individual missing from either source becomes a 0 rather
    than a dropped row. That is what makes the arm the whole condition.
    """
    index = pd.MultiIndex.from_product(
            [sorted(genes), sorted(individuals)],
            names=["ref_id", "individualID"])
    matrix = pd.DataFrame(index=index).reset_index()
    matrix = matrix.merge(tpm, on=["ref_id", "individualID"], how="left")
    matrix = matrix.merge(varcount, on=["ref_id", "individualID"], how="left")
    matrix["tpm"] = matrix["tpm"].fillna(0.0)
    matrix["varcount"] = matrix["varcount"].fillna(0).astype("int64")
    matrix["condition"] = matrix["individualID"].map(condition_of)
    return matrix


#----------------------------------------------------
# the tests
#----------------------------------------------------

def correlate(matrix, case, reference):
    """
    One row per gene: a Spearman per condition, then Fisher z on the
    difference.

    The Fisher z guards: both rho defined, both arms above
    MIN_INDIVIDUALS - 1, and both |rho| < 1, the last because arctanh(1) is
    infinite. The standard error is sqrt(1/(n1-3) + 1/(n2-3)), the Pearson
    form, without the 1.06 factor sometimes applied to Spearman.
    """
    records = []
    for ref_id, block in matrix.groupby("ref_id", sort=True):
        row = {"ref_id": ref_id}
        rho, sizes = {}, {}
        for label, prefix in ((case, "case"), (reference, "ref")):
            arm = block.loc[block["condition"] == label]
            testable = (len(arm) >= MIN_INDIVIDUALS
                        and arm["tpm"].nunique() > 1
                        and arm["varcount"].nunique() > 1)
            if testable:
                r, p = spearmanr(arm["tpm"], arm["varcount"])
            else:
                r, p = np.nan, np.nan
            rho[prefix], sizes[prefix] = r, len(arm)
            row[f"{prefix}_rho"] = r
            row[f"{prefix}_pvalue"] = p
            row[f"n_{prefix}"] = len(arm)

        r1, r2 = rho["case"], rho["ref"]
        n1, n2 = sizes["case"], sizes["ref"]
        if (pd.notna(r1) and pd.notna(r2)
                and n1 > 3 and n2 > 3
                and abs(r1) < 1 and abs(r2) < 1):
            z1, z2 = np.arctanh(r1), np.arctanh(r2)
            se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
            row["diff_pvalue"] = 2 * (1 - norm.cdf(abs((z1 - z2) / se)))
        else:
            row["diff_pvalue"] = np.nan
        records.append(row)

    tests = pd.DataFrame.from_records(records)
    if tests.empty:
        return tests

    # BH over the genes that produced a difference at all. A gene whose
    # correlation was not computable is absent from the ranking rather than
    # ranked as non-significant.
    tests["diff_adj_pvalue"] = np.nan
    testable = tests["diff_pvalue"].notna()
    if testable.any():
        tests.loc[testable, "diff_adj_pvalue"] = significance.benjamini_hochberg(
                tests.loc[testable, "diff_pvalue"].to_numpy())
    return tests


def summarise(tests, case, reference):
    """
    How many genes correlate in each condition, in both, and in neither.

    Nominal p-values, which is what the figure's fill colours encode.
    """
    case_sig = tests["case_pvalue"].notna() & (tests["case_pvalue"] < DIFF_ALPHA)
    ref_sig = tests["ref_pvalue"].notna() & (tests["ref_pvalue"] < DIFF_ALPHA)
    return {
        "both": int((case_sig & ref_sig).sum()),
        f"{case} only": int((case_sig & ~ref_sig).sum()),
        f"{reference} only": int((~case_sig & ref_sig).sum()),
        "neither": int((~case_sig & ~ref_sig).sum()),
        }


#----------------------------------------------------
# writing
#----------------------------------------------------

def check_reloadable(con, overwrite):
    counts = {table: con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
              for table in EXPRESSION_TABLES}
    if any(counts.values()) and not overwrite:
        raise SystemExit(
                f"tarva analysispipe expression-correlation: "
                f"{', '.join(f'{t} holds {n:,} row(s)' for t, n in counts.items())}. "
                f"Pass --overwrite to recompute.")


def prepare_tables(con):
    counts = {table: con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
              for table in EXPRESSION_TABLES}
    if any(counts.values()):
        for table in EXPRESSION_TABLES:
            con.execute(f"DROP TABLE IF EXISTS {table}")
        con.commit()
        schema.create_schema(con, EXPRESSION_TABLES)
        schema.report(f"--overwrite: dropped and recreated "
                      f"{', '.join(EXPRESSION_TABLES)}")


#----------------------------------------------------
# entry point
#----------------------------------------------------

def run(args):
    db_path = Path(args.db)
    if not db_path.is_file():
        raise SystemExit(
                f"tarva analysispipe expression-correlation: database not found: "
                f"{db_path}")

    con = schema.connect_for_build(db_path)
    try:
        schema.create_schema(con)
        check_reloadable(con, args.overwrite)

        case, reference, sizes = gene_tests.load_conditions(con, None)
        schema.report(f"conditions: case={case} (n={sizes[case]}), "
                      f"reference={reference} (n={sizes[reference]})")

        gene_set = load_gene_set(con)
        genes = set(gene_set["ref_id"])
        schema.report(f"genes: {len(genes):,} from gene_class")

        filtered = schema.get_meta(con, "gene_tests.subset_from_global_set")
        if filtered is not None and not int(filtered):
            schema.warn("expression-correlation: gene_class was built over ALL genes "
                        "because step 1a's global set was empty, so these correlations "
                        "describe the unfiltered gene set. run_meta records which set "
                        "was used.")

        table = schema.INDIVIDUAL_TABLE
        id_column, condition_column = schema.individual_columns(con)
        individuals = pd.read_sql_query(
                f'SELECT "{id_column}" AS individualID, '
                f'"{condition_column}" AS condition FROM {table}', con)
        condition_of = dict(zip(individuals["individualID"],
                                individuals["condition"]))

        tpm = load_tpm(con, genes)
        varcount = load_varcount(con, genes, args.min_depth)
        matrix = build_matrix(genes, set(individuals["individualID"]),
                              tpm, varcount, condition_of)
        schema.report(f"matrix: {len(matrix):,} (gene, individual) row(s); "
                      f"{len(tpm):,} with quantified expression, "
                      f"{len(varcount):,} with an eligible call")

        tests = correlate(matrix, case, reference).merge(
                gene_set, on="ref_id", how="left")

        n_case_rho = int(tests["case_rho"].notna().sum())
        n_ref_rho = int(tests["ref_rho"].notna().sum())
        schema.report(f"correlations: {n_case_rho:,} computable in {case}, "
                      f"{n_ref_rho:,} in {reference}, of {len(tests):,} gene(s)")

        if not n_case_rho or not n_ref_rho:
            # Every individual is present, so the arm size IS the cohort size:
            # a condition below MIN_INDIVIDUALS fails for every gene at once.
            short = [label for label in (case, reference)
                     if sizes[label] < MIN_INDIVIDUALS]
            if short:
                schema.warn(
                        f"expression-correlation: {', '.join(short)} has fewer than "
                        f"{MIN_INDIVIDUALS} individuals, so Spearman is not defined "
                        f"for ANY gene in that condition and no correlation can be "
                        f"computed. This is a property of the cohort -- the step is "
                        f"not applicable below {MIN_INDIVIDUALS} individuals per "
                        f"condition.")
            else:
                schema.warn(
                        f"expression-correlation: no gene had more than one distinct "
                        f"TPM and more than one distinct variant count in both "
                        f"conditions, so no correlation was computable. Check that "
                        f"load-expression ran and that --min-depth {args.min_depth} "
                        f"leaves calls in these genes.")

        n_diff = int(tests["diff_pvalue"].notna().sum())
        schema.report(f"between-condition tests: {n_diff:,} gene(s) with both "
                      f"correlations defined")

        if n_diff:
            comparable = tests.loc[tests["diff_pvalue"].notna()].rename(
                    columns={"diff_pvalue": significance.RAW,
                             "diff_adj_pvalue": significance.ADJUSTED})
            _, n_diff_sig = significance.choose_basis(
                    comparable, DIFF_ALPHA, significance.REPORTING,
                    "between-condition difference",
                    downstream="Figure 4's outlined points")
        else:
            n_diff_sig = 0

        groups = summarise(tests, case, reference)
        schema.report("correlated at nominal p < "
                      f"{DIFF_ALPHA}: " + ", ".join(
                          f"{label} {count:,}" for label, count in groups.items()))

        prepare_tables(con)
        tests.reindex(columns=[
            "ref_id", "gene_group", "novel", "n_case", "case_rho", "case_pvalue",
            "n_ref", "ref_rho", "ref_pvalue", "diff_pvalue",
            "diff_adj_pvalue"]).to_sql(
                "expression_test", con, if_exists="append", index=False,
                chunksize=gene_tests.WRITE_CHUNK)
        con.commit()
        schema.report(f"expression_test: {len(tests):,} row(s) written")

        schema.create_indexes(con, EXPRESSION_TABLES)
        schema.set_meta(con, "expression_correlation", {
            "expression_correlation.case_condition": case,
            "expression_correlation.reference_condition": reference,
            "expression_correlation.min_depth": args.min_depth,
            "expression_correlation.min_individuals": MIN_INDIVIDUALS,
            "expression_correlation.diff_alpha": DIFF_ALPHA,
            "expression_correlation.n_genes": len(tests),
            "expression_correlation.n_case_rho": n_case_rho,
            "expression_correlation.n_ref_rho": n_ref_rho,
            "expression_correlation.n_diff_tested": n_diff,
            "expression_correlation.n_diff_significant": n_diff_sig,
            "expression_correlation.group_counts": ",".join(
                    f"{k}={v}" for k, v in groups.items()),
            "expression_correlation.output_dir": (str(Path(args.output).resolve())
                                                  if args.output else ""),
            })
    finally:
        con.close()

    if args.output:
        directory = Path(args.output)
        directory.mkdir(parents=True, exist_ok=True)
        # The condition labels ride along as columns so the figure can label
        # its axes without being told which group was the case arm.
        out = tests.copy()
        out["case_condition"] = case
        out["ref_condition"] = reference
        gene_tests.write_csv(out, directory, RELATIONSHIP_CSV)

    schema.report("")
    schema.report("expression-correlation complete")
    schema.report(f"  genes                     {len(tests):,}")
    schema.report(f"  rho computable in {case:<9} {n_case_rho:,}")
    schema.report(f"  rho computable in {reference:<9} {n_ref_rho:,}")
    schema.report(f"  both conditions           {n_diff:,}")
    schema.report(f"  difference significant    {n_diff_sig:,}")
    for label, count in groups.items():
        schema.report(f"  {label:<25} {count:,}")
    schema.report(f"Wrote {db_path}")
    return 0
