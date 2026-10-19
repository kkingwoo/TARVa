"""
analysispipe gene-tests: the per-gene and per-call-type tests at the global
level, and the condition-enriched gene subsets.

Four things in one command, since each reads what the one before it wrote.
Every intermediate quantity is stored beside its verdict rather than
discarded, so each number is checkable without a re-run.

    `tarva analysispipe gene-tests --db <TARVa.db> [--case-condition <label>]
        [--min-depth 20] [--group-threshold 0.75] [--output <dir>] [--overwrite]`

Input:  a database with annotate, load-samples and load-calls done; it reads
        `call`, `site`, `gene` and `individual`. --case-condition names the
        case group, and without it the alphabetically first label is used.
Output: the gene_test, gene_calltype, gene_class and calltype_test tables;
        calltype_tests.csv and calltype_counts.csv (Figure 3) when --output
        is given; and gene_tests.* keys in run_meta.

In order:

  **Entry criterion.** A gene is analysed only if there is at least one
      condition in which EVERY individual carries at least one eligible call
      (dp >= --min-depth) in that gene -- any call type. The unit is the
      individual, not the call: one individual holding every call in a gene
      says nothing about whether a within-gene comparison has anything to
      compare. Genes where neither condition is complete are dropped before
      anything else, which is why one arm of every test is fully observed.

  **1. Within genes, ALL CALL TYPES POOLED.** One test per GENE. Per
      individual, per site, proportion = ad/dp; summed per individual over
      every site in the gene AND every call type, into one number each.
      Individuals of either condition with no eligible call enter as 0
      (n_*_observed and n_*_padded record how many of each, so a result
      resting mostly on imputed zeros is visible at a glance). Two-sided
      Mann-Whitney U, then rbc = 1 - 2U / (n_case * n_ref), then
      Benjamini-Hochberg across genes.

      No call-type dimension: the question is which GENES differ.
      gene_calltype carries the per-(gene, call type) aggregates the local
      level needs, with no test attached.

      The GLOBAL GENE SET is then taken ADAPTIVELY: if more than half the
      starting gene set clears adjusted p < 0.05, those genes are the set;
      otherwise the step says so and falls back to raw p < 0.05.
      run_meta's gene_tests.global_set_basis records which one fired.

  **2. Across genes, PER CALL TYPE.** The same test asked of a call TYPE
      rather than a gene. Per individual, every length-weighted proportion
      (ad/dp) x pw they carry of that call type -- across every position and
      every gene in the set.
      pw = 1/gene_span.

      Over the SAME completeness-filtered genes, independent of which genes
      the per-gene test selected: "across all genes" means every gene that
      survived filtering.

  **3. Condition-enriched subsets.** For global-set genes only:
      alternate-allele counts summed per condition over the gene's total.
      75-100% in one condition assigns that condition's enriched subset,
      25-75% assigns Shared. gene_class also carries edit_dist, the
      majority-rule label the site-uniqueness step groups by.

  **4. Call-type load per individual.** Over the same gene set: the number of
      alternate-allele SITES each individual carries of each call type. A
      count, not a proportion, since step 2 already asked whether the
      proportions differ. One row is one point on Figure 3's boxplot, so it
      goes to CSV rather than to a table -- there is no verdict attached, and
      `COUNT(*)` against call, site and gene_class reproduces it at any time.

The global gene set needs no table of its own: it is
`SELECT ref_id FROM gene_test WHERE p_value < 0.05`, one row per gene, and
gene_class holds a row only for those genes.

"""

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu

from tarva.analysispipe import significance
from tarva.dbpipe import schema

GENE_TEST_TABLES = ("gene_test", "gene_calltype", "gene_class", "calltype_test")

# calltype_test.gene_set: which gene set a per-call-type row describes.
#
# "complete" -- every gene passing the entry criterion, which is what "across
# all genes" means: all genes that survived filtering, not the subset the
# per-gene test then selected. The two are independent questions asked of one
# filtered gene set: 1.) which genes differ between two conditions for all call types, and which
# call types differ between conditions across all genes. Both provide a separate set of information
# and neither depends on the other.
#
# Figure 3 is the result of the per-call-type test, which is
# a global characterisation.
GENE_SET_COMPLETE = "complete"

CALLTYPE_TESTS_CSV = "calltype_tests.csv"
CALLTYPE_COUNTS_CSV = "calltype_counts.csv"

DEFAULT_MIN_DEPTH = 20
DEFAULT_GROUP_THRESHOLD = 0.75
REQUIRED_CONDITIONS = 2

# The threshold a within-gene test must clear for its gene to enter the global
# gene set. Not a flag: one value, and unlike the 75% subset cutoff there is
# no sweep over alternatives.
#
# Which p-value it applies to is decided by significance.choose_basis: this is
# a GATING step, so it falls back to raw p when too little survives BH.
GLOBAL_ALPHA = 0.05
WRITE_CHUNK = 100000


def add_arguments(parser):
    parser.add_argument("--db", required=True,
                        help="Database with annotate, load-samples and load-calls done")
    parser.add_argument("--case-condition", default=None, metavar="LABEL",
                        help="Which condition label is the case group, filling the "
                             "n_case_* columns (default: the alphabetically first of the "
                             "two). The choice is recorded in run_meta and only affects "
                             "labelling and the sign of rbc, not the p-value")
    parser.add_argument("--min-depth", type=int, default=DEFAULT_MIN_DEPTH,
                        help=f"A call counts toward completeness and toward the test only "
                             f"at or above this depth (default: {DEFAULT_MIN_DEPTH})")
    parser.add_argument("--group-threshold", type=float, default=DEFAULT_GROUP_THRESHOLD,
                        help=f"Share of a gene's total at which a condition is called "
                             f"enriched (default: {DEFAULT_GROUP_THRESHOLD})")
    parser.add_argument("--output", default=None, metavar="DIR",
                        help="Directory for the CSVs the figures read: the per-call-type "
                             "tests and each individual's call-type site counts. "
                             "Created if absent. Without it, everything still goes to "
                             "the database and the counts are still reported")
    parser.add_argument("--overwrite", action="store_true",
                        help="Drop and recompute gene_test, gene_class and calltype_test "
                             "instead of failing")


#----------------------------------------------------
# inputs
#----------------------------------------------------

def load_conditions(con, case_condition):
    """
    (case_label, ref_label, {label: n_individuals}). Exactly two are required.

    The condition column's name comes from run_meta, since `individual`
    carries the metadata CSV's own headers, and is quoted because a CSV may
    legitimately have a column called `group`.
    """
    table = schema.INDIVIDUAL_TABLE
    schema.require_table(con, table, "tarva dbpipe load-samples")
    _, condition_column = schema.individual_columns(con)
    rows = con.execute(
            f'SELECT "{condition_column}", COUNT(*) FROM {table} '
            f'GROUP BY "{condition_column}" ORDER BY "{condition_column}"'
            ).fetchall()
    sizes = {str(label): int(count) for label, count in rows}
    if len(sizes) != REQUIRED_CONDITIONS:
        raise SystemExit(
                f"tarva analysispipe gene-tests: `individual` holds "
                f"{len(sizes)} condition(s) "
                f"({', '.join(sizes)}), and exactly {REQUIRED_CONDITIONS} are required. "
                f"Every test here is a two-sample comparison.")

    labels = sorted(sizes)
    if case_condition is None:
        case = labels[0]
    elif case_condition in sizes:
        case = case_condition
    else:
        raise SystemExit(
                f"tarva dbpipe stats: --case-condition {case_condition!r} is not a "
                f"condition in `individual`. Found: {', '.join(labels)}")
    reference = [label for label in labels if label != case][0]
    return case, reference, sizes


def load_calls(con, min_depth):
    """
    One row per (ref_id, a1type, individualID, site) with its a1prop share.

    Joined here rather than aggregated in SQL because the completeness filter,
    the padding and the test all need the per-sample values, and the eligible
    set is small once --min-depth and a non-NULL a1type are applied.
    """
    schema.require_table(con, "call", "tarva dbpipe load-calls")
    schema.require_table(con, "gene", "tarva dbpipe annotate")
    table = schema.INDIVIDUAL_TABLE
    id_column, condition_column = schema.individual_columns(con)
    query = f"""
        SELECT s.ref_id, c.a1type, c.individualID, c.a1_ad, c.dp, g.gene_span,
               p."{condition_column}" AS condition
        FROM call AS c
        JOIN site AS s ON s.chrom = c.chrom AND s.pos = c.pos
        JOIN gene AS g ON g.ref_id = s.ref_id
        JOIN {table} AS p ON p."{id_column}" = c.individualID
        WHERE s.ref_id IS NOT NULL
          AND c.a1type IS NOT NULL
          AND c.dp >= ?
    """
    frame = pd.read_sql_query(query, con, params=(int(min_depth),))
    if frame.empty:
        raise SystemExit(
                f"tarva dbpipe stats: no call has dp >= {min_depth}, a non-NULL a1type "
                f"and a gene. Nothing to test.")
    # Two quantities, used at different scales. WITHIN a gene (the per-gene
    # test) the proportion is plain ad/dp -- no length term, since every site
    # in the gene would carry the same one. ACROSS genes (the per-call-type
    # test) each proportion is multiplied by a position weight pw = 1/length,
    # which removes the gene-length bias.
    #
    # pw = 1/gene_span, where gene_span is the transcript length: the distance
    # from the gene's earliest transcript start to its latest transcript end,
    # introns included, counted in BASES (gene_end - gene_start + 1). Every
    # call in a gene shares the weight.
    frame["proportion"] = frame["a1_ad"] / frame["dp"]
    frame["weighted"] = frame["proportion"] / frame["gene_span"]
    return frame


#----------------------------------------------------
# step 1: completeness
#----------------------------------------------------

def complete_genes(frame, sizes):
    """
    {ref_id: complete_in} for genes that pass the entry criterion.

    A gene qualifies when at least one condition has EVERY ONE of its
    individuals carrying at least one eligible call in that gene.

    The unit is the INDIVIDUAL, not the call.

    Where both conditions are complete, both labels are recorded.
    """
    per_condition = (frame.groupby(["ref_id", "condition"])["individualID"]
                     .nunique().rename("n_individuals").reset_index())
    per_condition["full"] = per_condition["condition"].map(sizes)
    complete = per_condition.loc[
            per_condition["n_individuals"] >= per_condition["full"]]
    qualifying = {}
    for ref_id, condition in zip(complete["ref_id"], complete["condition"]):
        qualifying.setdefault(ref_id, []).append(condition)
    return {ref_id: ",".join(sorted(labels))
            for ref_id, labels in qualifying.items()}


#----------------------------------------------------
# 1. the within-gene test, all call types pooled
#----------------------------------------------------

def run_tests(frame, qualifying, case, reference, sizes):
    """
    One row per qualifying GENE: the padded Mann-Whitney and rbc, with every
    call type pooled.

    ONE OBSERVATION IS ONE INDIVIDUAL. All of that individual's unweighted
    ad/dp proportions within the gene -- across every position AND EVERY CALL
    TYPE -- sum into a single number, and the arms are those per-individual
    sums, zero-padded to cohort size. An individual with no eligible call in
    the gene enters as 0; n_*_observed and n_*_padded record how many of each,
    so a result resting mostly on imputed zeros is visible at a glance.

    n_sites_* are raw SITE counts per condition, carried beside the sums
    because they describe how much data the gene held, which a p-value does
    not.

    p_value is NaN when both arms are all-zero, which is a different statement
    from a p-value of 1 and is stored as NULL.
    """
    per_individual = (frame.groupby(["ref_id", "condition", "individualID"])
                      ["proportion"].sum().rename("total").reset_index())

    sites = frame.groupby(["ref_id", "condition"]).size().rename("n").reset_index()
    site_counts = {(ref_id, condition): int(n) for ref_id, condition, n
                   in zip(sites["ref_id"], sites["condition"], sites["n"])}

    records = []
    for ref_id, block in per_individual.groupby("ref_id", sort=True):
        complete_in = qualifying.get(ref_id)
        if complete_in is None:
            continue
        arms, counts = {}, {}
        for label in (case, reference):
            observed = block.loc[block["condition"] == label, "total"].to_numpy()
            padding = max(sizes[label] - observed.size, 0)
            arms[label] = np.concatenate([observed, np.zeros(padding)])
            counts[label] = (observed.size, padding)

        case_arm, ref_arm = arms[case], arms[reference]
        if case_arm.size and ref_arm.size and (case_arm.any() or ref_arm.any()):
            u_stat, p_value = mannwhitneyu(case_arm, ref_arm, alternative="two-sided")
            rbc = 1.0 - (2.0 * u_stat) / (case_arm.size * ref_arm.size)
        else:
            u_stat, p_value, rbc = np.nan, np.nan, np.nan

        records.append({
            "ref_id": ref_id, "complete_in": complete_in,
            "u_stat": u_stat, "p_value": p_value, "rbc": rbc,
            "n_case_observed": counts[case][0], "n_case_padded": counts[case][1],
            "n_ref_observed": counts[reference][0], "n_ref_padded": counts[reference][1],
            "n_sites_case": site_counts.get((ref_id, case), 0),
            "n_sites_ref": site_counts.get((ref_id, reference), 0),
            })

    tests = pd.DataFrame.from_records(records)
    if tests.empty:
        return tests
    tests["adj_pvalue"] = significance.benjamini_hochberg(
            tests["p_value"].to_numpy())
    return tests


def aggregate_calltypes(frame, case, reference):
    """
    One row per (ref_id, a1type): the per-condition aggregates local stage 4
    needs, and no test.

    Stage 4 runs a two-sample t-test on per-individual Sigma ad / Sigma dp
    within a gene.
    """
    stats = frame.groupby(["ref_id", "a1type", "condition"]).agg(
            sum_ad=("a1_ad", "sum"),
            sum_dp=("dp", "sum"),
            n_sites=("dp", "size"),
            n_indiv=("individualID", "nunique"),
            ).reset_index()

    lookup = {(row.ref_id, row.a1type, row.condition):
              (int(row.sum_ad), int(row.sum_dp), int(row.n_sites), int(row.n_indiv))
              for row in stats.itertuples()}

    records = []
    for ref_id, a1type in sorted({(row.ref_id, row.a1type)
                                  for row in stats.itertuples()}):
        case_ad, case_dp, case_sites, case_indiv = lookup.get(
                (ref_id, a1type, case), (0, 0, 0, 0))
        ref_ad, ref_dp, ref_sites, ref_indiv = lookup.get(
                (ref_id, a1type, reference), (0, 0, 0, 0))
        records.append({
            "ref_id": ref_id, "a1type": a1type,
            "n_sites_case": case_sites, "n_sites_ref": ref_sites,
            "n_indiv_case": case_indiv, "n_indiv_ref": ref_indiv,
            "sum_ad_case": case_ad, "sum_dp_case": case_dp,
            "sum_ad_ref": ref_ad, "sum_dp_ref": ref_dp,
            })
    return pd.DataFrame.from_records(records)


#----------------------------------------------------
# 2. the across-gene test, per call type
#----------------------------------------------------

def run_calltype_tests(frame, case, reference, gene_set, sizes):
    """
    One row per call type: the across-gene Mann-Whitney on length-weighted
    proportions, pooled over `gene_set`.

    The within-gene test asks whether one GENE's proportions differ between
    conditions. This asks whether a CALL TYPE does.

    ONE OBSERVATION IS ONE INDIVIDUAL. Every weighted proportion that
    individual carries of that call type -- across every position AND every
    gene in the set -- is summed into a single number for them, and the two
    arms are those per-individual sums. So the arms are the size of the
    cohort, not of the call set, and they are zero-padded as the within-gene
    test's are: an individual with no call of that type anywhere in the set
    carries none of it, which is a measurement and not a gap.

    pw = 1/gene_span is what lets a site in a 4 kb gene and a site in a 430 kb
    gene contribute comparably to that sum, instead of the long gene setting
    the total on its own.

    p_value is NaN when both arms are all-zero, and is stored as NULL.
    """
    per_individual = (frame.groupby(["a1type", "condition", "individualID"])
                      ["weighted"].sum().rename("total").reset_index())

    records = []
    for a1type, block in per_individual.groupby("a1type", sort=True):
        arms, counts = {}, {}
        for label in (case, reference):
            observed = block.loc[block["condition"] == label, "total"].to_numpy()
            padding = max(sizes[label] - observed.size, 0)
            arms[label] = np.concatenate([observed, np.zeros(padding)])
            counts[label] = (observed.size, padding)

        case_arm, ref_arm = arms[case], arms[reference]
        if case_arm.size and ref_arm.size and (case_arm.any() or ref_arm.any()):
            u_stat, p_value = mannwhitneyu(case_arm, ref_arm, alternative="two-sided")
            rbc = 1.0 - (2.0 * u_stat) / (case_arm.size * ref_arm.size)
        else:
            u_stat, p_value, rbc = np.nan, np.nan, np.nan

        records.append({
            "gene_set": gene_set, "a1type": a1type,
            "u_stat": u_stat, "p_value": p_value, "rbc": rbc,
            "n_case_observed": counts[case][0], "n_case_padded": counts[case][1],
            "n_ref_observed": counts[reference][0],
            "n_ref_padded": counts[reference][1],
            })

    tests = pd.DataFrame.from_records(records)
    if tests.empty:
        return tests
    tests["adj_pvalue"] = significance.benjamini_hochberg(
            tests["p_value"].to_numpy())
    return tests


#----------------------------------------------------
# 3. condition-enriched gene subsets
#----------------------------------------------------

def assign_groups(n_case, n_ref, case, reference, threshold):
    """
    Each gene's group label at `threshold`, given alternate-allele counts.

    """
    n_case = np.asarray(n_case, dtype=float)
    n_ref = np.asarray(n_ref, dtype=float)
    total = n_case + n_ref
    present = total > 0
    case_share = np.divide(n_case, total, out=np.zeros_like(total), where=present)
    ref_share = np.divide(n_ref, total, out=np.zeros_like(total), where=present)
    return np.where(
            case_share >= threshold, f"{case}_enriched",
            np.where(ref_share >= threshold, f"{reference}_enriched", "Shared"))


def classify_genes(frame, case, reference, threshold):
    """
    Over the GLOBAL GENE SET only: the number of alternate alleles per
    individual summed across each condition, divided by the gene's total.
    75-100% in one condition assigns that condition's enriched subset;
    25-75% assigns Shared.

    Counts of alternate alleles, not summed allele depths.

    edit_dist is the second, looser classification of the same gene -- the
    per-gene distribution label the site-uniqueness step groups by. Whichever
    condition holds MORE THAN HALF gives that condition's *_enriched label,
    and an exact tie gives Both. Three categories, produced as three.

    """
    columns = ["ref_id", "condition"]
    if frame.empty:
        return pd.DataFrame(columns=["ref_id", "edit_dist", "gene_group",
                                     "n_calls_case", "n_calls_ref"])

    calls = frame.groupby(columns).size().rename("n_calls").reset_index()
    wide = calls.pivot(index="ref_id", columns="condition",
                       values="n_calls").fillna(0)
    for label in (case, reference):
        if label not in wide.columns:
            wide[label] = 0

    out = pd.DataFrame(index=wide.index)
    out["n_calls_case"] = wide[case].astype("int64")
    out["n_calls_ref"] = wide[reference].astype("int64")

    out["gene_group"] = assign_groups(
            out["n_calls_case"], out["n_calls_ref"], case, reference, threshold)

    # edit_dist needs the shares themselves, which assign_groups does not
    # return; it is a per-gene distribution label, not a subset assignment.
    total = out["n_calls_case"] + out["n_calls_ref"]
    case_share = np.where(total > 0, out["n_calls_case"] / total.where(total > 0, 1), 0.0)
    ref_share = 1.0 - case_share

    out["edit_dist"] = np.select(
            [case_share > 0.5, ref_share > 0.5],
            [f"{case}_enriched", f"{reference}_enriched"],
            default="Both")
    return out.reset_index()


#----------------------------------------------------
# 4. call-type load per individual in the global set
#----------------------------------------------------

def count_calltype_sites(frame, genes):
    """
    One row per (call type, individual): how many alternate-allele sites that
    individual carries of that type, over `genes`.

    """
    subset = frame.loc[frame["ref_id"].isin(genes)]
    counts = (subset.groupby(["a1type", "individualID", "condition"])
              .size().rename("total sites").reset_index())
    counts = counts.rename(columns={"a1type": "variant_call"})
    return counts.reindex(columns=["variant_call", "total sites", "condition",
                                   "individualID"])


#----------------------------------------------------
# writing
#----------------------------------------------------

def write_csv(frame, directory, name):
    path = directory / name
    frame.to_csv(path, index=False)
    schema.report(f"wrote {path} ({len(frame):,} row(s))")
    return path


def check_reloadable(con, overwrite):
    counts = {table: con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
              for table in GENE_TEST_TABLES}
    if any(counts.values()) and not overwrite:
        raise SystemExit(
                f"tarva analysispipe gene-tests: "
                f"{', '.join(f'{t} holds {n:,} row(s)' for t, n in counts.items())}. "
                f"Pass --overwrite to recompute them.")


def prepare_tables(con):
    counts = {table: con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
              for table in GENE_TEST_TABLES}
    if any(counts.values()):
        for table in GENE_TEST_TABLES:
            con.execute(f"DROP TABLE IF EXISTS {table}")
        con.commit()
        schema.create_schema(con, GENE_TEST_TABLES)
        schema.report(f"--overwrite: dropped and recreated "
                      f"{', '.join(GENE_TEST_TABLES)}")


#----------------------------------------------------
# entry point
#----------------------------------------------------

def run(args):
    db_path = Path(args.db)
    if not db_path.is_file():
        raise SystemExit(
                f"tarva analysispipe gene-tests: database not found: {db_path}")

    con = schema.connect_for_build(db_path)
    try:
        schema.create_schema(con)
        check_reloadable(con, args.overwrite)

        case, reference, sizes = load_conditions(con, args.case_condition)
        schema.report(f"conditions: case={case} (n={sizes[case]}), "
                      f"reference={reference} (n={sizes[reference]})")

        frame = load_calls(con, args.min_depth)
        schema.report(f"eligible calls: {len(frame):,} at dp >= {args.min_depth}, "
                      f"over {frame['ref_id'].nunique():,} gene(s)")

        genes_before = frame["ref_id"].nunique()
        qualifying = complete_genes(frame, sizes)
        schema.report(f"completeness: {len(qualifying):,} of {genes_before:,} gene(s) "
                      f"have a condition in which every individual carries a call")
        if not qualifying:
            raise SystemExit(
                    f"tarva analysispipe gene-tests: no gene has even one condition in "
                    f"which EVERY individual carries an eligible call, so nothing is "
                    f"testable. With {sizes[case]} and {sizes[reference]} individuals "
                    f"that criterion is strict; check the eligible-call count above.")

        # Completeness decides which GENES are tested.
        analysed = frame.loc[frame["ref_id"].isin(set(qualifying))]
        schema.report(f"analysed: {len(analysed):,} of {len(frame):,} eligible call(s) "
                      f"in the surviving gene(s)")

        tests = run_tests(analysed, qualifying, case, reference, sizes)
        testable = int(tests["p_value"].notna().sum())
        schema.report(f"per-gene tests: {len(tests):,} with all call types "
                      f"pooled, {testable:,} with a defined p-value")

        # No test here -- these are the local level's inputs, at its unit.
        calltype_rows = aggregate_calltypes(analysed, case, reference)
        schema.report(f"gene_calltype: {len(calltype_rows):,} (gene, call type) "
                      f"aggregate(s)")

        # Retention: only genes clearing GLOBAL_ALPHA in their within-gene
        # test enter the GLOBAL GENE SET, and the subset assignment onward
        # operates on that set -- not on every gene that passed completeness.
        # The set itself needs no table: it is
        # `SELECT DISTINCT ref_id FROM gene_test WHERE p_value < 0.05`.
        basis, _ = significance.choose_basis(
                tests, GLOBAL_ALPHA, significance.GATING, "the per-gene test",
                downstream="the condition-enriched subset assignment")
        # basis is None when neither adjusted nor nominal p reached alpha in a
        # majority of tests, which selects no genes rather than keeping the
        # few that happened to clear a nominal threshold.
        global_genes = (set() if basis is None else
                        set(tests.loc[tests[basis] < GLOBAL_ALPHA, "ref_id"]))
        schema.report(f"global gene set: {len(global_genes):,} gene(s) "
                      f"({basis or 'no qualifying basis'})")

        # The per-call-type test runs on the completeness-filtered set, NOT on
        # the per-gene test's selection: across all genes means every gene
        # that survived filtering. It therefore needs no gene set resolved for
        # it and cannot be affected by the per-gene test selecting nothing.
        calltypes = run_calltype_tests(analysed, case, reference,
                                       GENE_SET_COMPLETE, sizes)
        schema.report(f"per-call-type tests: {len(calltypes):,} call type(s) over "
                      f"{len(analysed):,} weighted proportion(s) summed into "
                      f"{sizes[case]} vs {sizes[reference]} per-individual total(s), "
                      f"{GENE_SET_COMPLETE} gene set")
        calltype_basis, n_calltype_sig = significance.choose_basis(
                calltypes, GLOBAL_ALPHA, significance.REPORTING,
                "the per-call-type test",
                downstream="Figure 3's significance annotation")
        # Stored so Figure 3 marks the column this step actually used, rather
        # than assuming adjusted. Empty when neither basis qualified.
        calltypes["basis"] = calltype_basis or ""

        # One resolution for both the subsets and the call-type load, so they
        # describe the same genes.
        subset_genes, filtered = significance.resolve_gene_set(
                global_genes, set(analysed["ref_id"]),
                "the subsets and the call-type load", "the per-gene test")

        classes = classify_genes(frame.loc[frame["ref_id"].isin(subset_genes)],
                                 case, reference, args.group_threshold)
        schema.report(f"gene groups: {len(classes):,} gene(s) assigned at "
                      f"{args.group_threshold:.0%}")

        counts = count_calltype_sites(frame, subset_genes)
        n_subset_calls = int(counts["total sites"].sum())
        schema.report(f"call-type load: {n_subset_calls:,} call(s) over "
                      f"{len(subset_genes):,} gene(s), {len(counts):,} "
                      f"(call type, individual) row(s)")

        prepare_tables(con)
        tests.reindex(columns=[
            "ref_id", "complete_in", "u_stat", "p_value", "adj_pvalue", "rbc",
            "n_case_observed", "n_case_padded", "n_ref_observed", "n_ref_padded",
            "n_sites_case", "n_sites_ref"]).to_sql(
                "gene_test", con, if_exists="append", index=False, chunksize=WRITE_CHUNK)
        con.commit()
        schema.report(f"gene_test: {len(tests):,} row(s) written")

        calltype_rows.reindex(columns=[
            "ref_id", "a1type", "n_sites_case", "n_sites_ref", "n_indiv_case",
            "n_indiv_ref", "sum_ad_case", "sum_dp_case", "sum_ad_ref",
            "sum_dp_ref"]).to_sql(
                "gene_calltype", con, if_exists="append", index=False,
                chunksize=WRITE_CHUNK)
        con.commit()
        schema.report(f"gene_calltype: {len(calltype_rows):,} row(s) written")

        classes.reindex(columns=[
            "ref_id", "edit_dist", "gene_group", "n_calls_case",
            "n_calls_ref"]).to_sql(
                "gene_class", con, if_exists="append", index=False, chunksize=WRITE_CHUNK)
        con.commit()
        schema.report(f"gene_class: {len(classes):,} row(s) written")

        calltypes.reindex(columns=[
            "gene_set", "a1type", "u_stat", "p_value", "adj_pvalue", "basis",
            "rbc", "n_case_observed", "n_case_padded", "n_ref_observed",
            "n_ref_padded"]).to_sql(
                "calltype_test", con, if_exists="append", index=False,
                chunksize=WRITE_CHUNK)
        con.commit()
        schema.report(f"calltype_test: {len(calltypes):,} row(s) written")

        schema.create_indexes(con, GENE_TEST_TABLES)

        if args.output:
            directory = Path(args.output)
            directory.mkdir(parents=True, exist_ok=True)
            write_csv(calltypes, directory, CALLTYPE_TESTS_CSV)
            write_csv(counts, directory, CALLTYPE_COUNTS_CSV)

        group_counts = classes["gene_group"].value_counts()
        schema.set_meta(con, "gene_tests", {
            "gene_tests.case_condition": case,
            "gene_tests.reference_condition": reference,
            "gene_tests.n_case": sizes[case],
            "gene_tests.n_reference": sizes[reference],
            "gene_tests.min_depth": args.min_depth,
            "gene_tests.global_alpha": GLOBAL_ALPHA,
            "gene_tests.group_threshold": args.group_threshold,
            "gene_tests.n_eligible_calls": len(frame),
            "gene_tests.n_genes_considered": genes_before,
            "gene_tests.n_genes_complete": len(qualifying),
            "gene_tests.n_tests": len(tests),
            "gene_tests.n_tests_defined": testable,
            "gene_tests.n_calltype_rows": len(calltype_rows),
            "gene_tests.n_genes_global_set": len(global_genes),
            "gene_tests.global_set_basis": basis or "none",
            "gene_tests.calltype_basis": calltype_basis or "none",
            "gene_tests.calltype_gene_set": GENE_SET_COMPLETE,
            "gene_tests.n_calltypes_tested": len(calltypes),
            "gene_tests.n_calltypes_adjusted_sig": n_calltype_sig,
            # 0 means the subset assignment fell back to every available gene
            # because the global set was empty, so gene_class is not a
            # filtered result.
            "gene_tests.subset_from_global_set": int(filtered),
            "gene_tests.n_genes_in_subsets": len(subset_genes),
            "gene_tests.n_calls_in_subsets": n_subset_calls,
            "gene_tests.output_dir": (str(Path(args.output).resolve())
                                      if args.output else ""),
            "gene_tests.n_tests_adjusted_sig": int(
                    (tests["adj_pvalue"] < GLOBAL_ALPHA).sum()),
            "gene_tests.n_genes_classified": len(classes),
            "gene_tests.group_sizes": ",".join(
                    f"{k}={v}" for k, v in group_counts.items()),
            })
    finally:
        con.close()

    schema.report("")
    schema.report("gene-tests complete")
    schema.report(f"  genes complete            {len(qualifying):,} of {genes_before:,}")
    schema.report(f"  genes tested              {len(tests):,}")
    schema.report(f"  with a defined p-value    {testable:,}")
    schema.report(f"  (gene, call type) rows    {len(calltype_rows):,}")
    bases = {significance.ADJUSTED: "adjusted", significance.RAW: "nominal"}
    schema.report(f"  global gene set           {len(global_genes):,} "
                  f"({bases.get(basis, 'no qualifying basis')}"
                  f"{'' if basis is None else f' p < {GLOBAL_ALPHA}'})")
    # raw p < 0.05 first: it is the threshold the global set falls back to
    for column, level, label in (("p_value", GLOBAL_ALPHA, "raw p < 0.05"),
                                 ("adj_pvalue", GLOBAL_ALPHA, "adj p < 0.05"),
                                 ("p_value", 0.01, "raw p < 0.01"),
                                 ("p_value", 0.001, "raw p < 0.001"),
                                 ("p_value", 0.2, "raw p < 0.2")):
        schema.report(f"  {label:<25} "
                      f"{int((tests[column] < level).sum()):,}")
    schema.report(f"  smallest raw p            "
                  f"{significance.format_p(tests['p_value'].min())}")
    schema.report("")
    schema.report(f"  per call type: {len(calltypes):,} call type(s), {n_calltype_sig:,} "
                  f"significant on {bases.get(calltype_basis, 'no qualifying basis')}"
                  f"{'' if calltype_basis is None else f' p < {GLOBAL_ALPHA}'}")
    for row in calltypes.sort_values("p_value").itertuples():
        schema.report(f"    {row.a1type:<8} "
                      f"nominal p {significance.format_p(row.p_value):>10}   "
                      f"adj p {significance.format_p(row.adj_pvalue):>10}   "
                      f"rbc {row.rbc:>7.3f}   "
                      f"observed {row.n_case_observed}/{row.n_ref_observed}")
    schema.report("")
    if not filtered:
        schema.report(f"  subsets computed over    ALL {len(subset_genes):,} genes "
                      f"(global set was empty)")
    for group, count in group_counts.items():
        schema.report(f"  {group:<25} {count:,}")
    schema.report(f"  call-type load            {n_subset_calls:,} call(s) in "
                  f"{len(counts):,} (call type, individual) row(s)")
    schema.report(f"Wrote {db_path}")
    return 0
