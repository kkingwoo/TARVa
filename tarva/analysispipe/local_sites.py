"""
analysispipe local-sites: the local level, stages 4 to 6.

Having asked which genes and which call types differ, ask which (gene, call
type) pairs differ, and then which individual positions within them.

  STAGE 4  Per (gene, call type), the aggregated proportion per condition --
           sum(ad) over sum(dp) across every alternate-allele site and every
           carrier of that condition. Two-sample Welch t-test on the
           per-individual aggregates, then Benjamini-Hochberg.

           ARMS ARE CARRIERS, NOT THE COHORT. Steps 1a and 1b zero-pad an
           absent individual because their question is about the whole
           condition. Here the retention rule does that work instead: a pair
           is kept only when at least half the individuals of one condition
           carry the call type, so a surviving arm is populated by
           construction and padding would only add zeros to both sides.

           Retention: adjusted-or-nominal p below LOCAL_ALPHA (0.01) AND at
           least half the individuals of one condition carrying that call
           type in that gene. The unit there is the individual, not the site:
           one individual holding fifteen sites is not fifteen individuals'
           worth of evidence.

  STAGE 5  Distinct individuals per condition carrying the pair, and the local
           distribution label from them:

             specific      one condition only, and ALL of its individuals
             low-specific  one condition only, but not all of them
             double        both conditions, one carrier count at least 1.5x
                           the other
             differential  both conditions, neither reaching 1.5x

  STAGE 6  Mann-Whitney over the length-weighted per-site proportions, pooled
           across the gene's positions and split by condition, with the
           condition whose mean is larger recorded.

           THIS UNIT IS THE (INDIVIDUAL, POSITION) PAIR, not the individual --
           see the note below.

  FINAL    One condition's carrier count above 1 and at least half again the
           other's, OR a nominal stage-6 p below NOMINAL_ALPHA.

A NOTE ON STAGE 6'S UNIT. It is the (individual, position) pair, not the
individual: proportions accumulate per position across the gene and one
Mann-Whitney runs over those values, so an individual carrying several
positions contributes several observations. Stage 4 already asks whether the
per-individual aggregate differs, so a per-individual stage 6 would be the
same question asked non-parametrically; per-position instead asks whether the
DISTRIBUTION ACROSS POSITIONS differs. The arms are therefore not independent
observations, and this p-value is used only as one of two alternative
retention criteria, never as a cohort-level significance claim.

    `tarva analysispipe local-sites --db <TARVa.db> [--min-depth 20]
        [--output <dir>] [--overwrite]`

Input:  a database with `gene-tests` done, for the gene set, the subset
        labels and the condition labels.
Output: the local_calltype_test and local_site tables,
        local_calltype_tests.csv (Table 2) and local_sites.csv (Table 3) when
        --output is given, and local_sites.* keys in run_meta.
"""

import math

from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu, ttest_ind

from tarva.analysispipe import gene_tests
from tarva.analysispipe import significance
from tarva.dbpipe import schema

LOCAL_TABLES = ("local_calltype_test", "local_site")

CALLTYPE_CSV = "local_calltype_tests.csv"
SITES_CSV = "local_sites.csv"

# Stage 4's retention threshold. Not a flag.
LOCAL_ALPHA = 0.01

# The looser threshold for counting genes with at least one call-type
# difference, and stage 6's final rule.
NOMINAL_ALPHA = 0.05

# "At least half again the other", i.e. 50% more.
DOMINANCE = 1.5

DEFAULT_MIN_DEPTH = 20

SPECIFIC, LOW_SPECIFIC = "specific", "low-specific"
DOUBLE, DIFFERENTIAL = "double", "differential"


def add_arguments(parser):
    parser.add_argument("--db", required=True,
                        help="Database with gene-tests done")
    parser.add_argument("--min-depth", type=int, default=DEFAULT_MIN_DEPTH,
                        help=f"A call enters the local stages only at or above this "
                             f"depth (default: {DEFAULT_MIN_DEPTH}). Should match the "
                             f"--min-depth gene-tests used")
    parser.add_argument("--output", default=None, metavar="DIR",
                        help="Directory for local_calltype_tests.csv (Table 2) and "
                             "local_sites.csv (Table 3). Created if absent")
    parser.add_argument("--overwrite", action="store_true",
                        help="Drop and recompute the local tables instead of failing")


#----------------------------------------------------
# inputs
#----------------------------------------------------

def load_calls(con, min_depth):
    """
    One row per (gene, call type, individual, position): ad, dp, gene span,
    condition.

    Restricted to the gene set steps 2a and 2b used, via gene_class, so the
    local level narrows the global result rather than reopening it.
    """
    schema.require_table(con, "gene_class", "tarva analysispipe gene-tests")
    schema.require_table(con, "call", "tarva dbpipe load-calls")
    table = schema.INDIVIDUAL_TABLE
    id_column, condition_column = schema.individual_columns(con)
    query = f"""
        SELECT s.ref_id        AS ref_id,
               c.a1type        AS a1type,
               c.individualID  AS individualID,
               s.chrom         AS chrom,
               s.pos           AS pos,
               c.a1_ad         AS a1_ad,
               c.dp            AS dp,
               g.gene_span     AS gene_span,
               p."{condition_column}" AS condition
        FROM call AS c
        JOIN site AS s ON s.chrom = c.chrom AND s.pos = c.pos
        JOIN gene_class AS k ON k.ref_id = s.ref_id
        JOIN gene AS g ON g.ref_id = s.ref_id
        JOIN {table} AS p ON p."{id_column}" = c.individualID
        WHERE c.a1type IS NOT NULL
          AND c.dp >= ?
    """
    frame = pd.read_sql_query(query, con, params=(int(min_depth),))
    if frame.empty:
        raise SystemExit(
                f"tarva analysispipe local-sites: no call in the gene_class gene set "
                f"has a call type and dp >= {min_depth}. Nothing to test.")
    # The same two quantities steps 1a and 1b use, for the same reasons.
    frame["proportion"] = frame["a1_ad"] / frame["dp"]
    frame["weighted"] = frame["proportion"] / frame["gene_span"]
    return frame


#----------------------------------------------------
# stage 4
#----------------------------------------------------

def majority(size):
    """At least half of `size` individuals, rounded up."""
    return math.ceil(size / 2)


def stage4(frame, case, reference, sizes):
    """
    One row per (gene, call type): per-individual aggregated proportions,
    Welch's t-test, and the per-condition carrier counts.

    The aggregate is sum(ad) / sum(dp) -- summed first, divided once. Not the
    mean of per-site proportions, which would weight a shallow site equally
    with a deep one.
    """
    per_individual = (frame.groupby(
            ["ref_id", "a1type", "condition", "individualID"])[["a1_ad", "dp"]]
            .sum().reset_index())
    per_individual["proportion"] = (per_individual["a1_ad"]
                                    / per_individual["dp"])

    records = []
    for (ref_id, a1type), block in per_individual.groupby(["ref_id", "a1type"],
                                                          sort=True):
        arms, counts, aggregate = {}, {}, {}
        for label, prefix in ((case, "case"), (reference, "ref")):
            arm = block.loc[block["condition"] == label]
            arms[prefix] = arm["proportion"].to_numpy()
            counts[prefix] = len(arm)
            total_dp = arm["dp"].sum()
            aggregate[prefix] = (float(arm["a1_ad"].sum() / total_dp)
                                 if total_dp else np.nan)

        if counts["case"] > 1 and counts["ref"] > 1:
            t_stat, p_value = ttest_ind(arms["case"], arms["ref"],
                                        equal_var=False)
        else:
            t_stat, p_value = np.nan, np.nan

        records.append(dict(
            ref_id=ref_id, a1type=a1type,
            n_case=counts["case"], n_ref=counts["ref"],
            prop_case=aggregate["case"], prop_ref=aggregate["ref"],
            t_stat=t_stat, p_value=p_value))

    tests = pd.DataFrame.from_records(records)
    if tests.empty:
        return tests
    tests["adj_pvalue"] = significance.benjamini_hochberg(
            tests["p_value"].to_numpy())
    # Half of one condition's individuals, whichever condition reaches it.
    tests["has_majority"] = (
            (tests["n_case"] >= majority(sizes[case]))
            | (tests["n_ref"] >= majority(sizes[reference])))
    return tests


#----------------------------------------------------
# stage 5
#----------------------------------------------------

def classify_local(tests, case, reference, sizes):
    """
    The local distribution label per (gene, call type), from the carrier
    counts alone.

    `specific` requires every individual of the single present condition, so
    it is the strongest statement available: this call type appears in that
    condition and in all of it.
    """
    n_case, n_ref = tests["n_case"], tests["n_ref"]
    both = (n_case > 0) & (n_ref > 0)
    larger = np.maximum(n_case, n_ref)
    smaller = np.minimum(n_case, n_ref)
    complete = ((n_ref == 0) & (n_case == sizes[case])) | \
               ((n_case == 0) & (n_ref == sizes[reference]))
    return np.select(
            [~both & complete,
             ~both,
             both & (larger >= DOMINANCE * smaller)],
            [SPECIFIC, LOW_SPECIFIC, DOUBLE],
            default=DIFFERENTIAL)


#----------------------------------------------------
# stage 6
#----------------------------------------------------

def stage6(frame, pairs, case, reference):
    """
    Per (gene, call type), a Mann-Whitney over the length-weighted per-site
    proportions, split by condition, plus which condition's mean was larger.

    `pairs` is the set of (ref_id, a1type) stage 4 retained -- the whole point
    of the cascade is that this runs on those alone.
    """
    records = []
    if not pairs:
        return pd.DataFrame(columns=["ref_id", "a1type", "prop_u_stat",
                                     "prop_p_value", "higher_prop"])

    wanted = frame.set_index(["ref_id", "a1type"]).index.isin(pairs)
    for (ref_id, a1type), block in frame.loc[wanted].groupby(
            ["ref_id", "a1type"], sort=True):
        case_arm = block.loc[block["condition"] == case, "weighted"].to_numpy()
        ref_arm = block.loc[block["condition"] == reference, "weighted"].to_numpy()
        if case_arm.size and ref_arm.size:
            u_stat, p_value = mannwhitneyu(case_arm, ref_arm,
                                           alternative="two-sided")
            higher = case if case_arm.mean() > ref_arm.mean() else reference
        else:
            u_stat, p_value = np.nan, np.nan
            higher = case if case_arm.size else reference
        records.append(dict(ref_id=ref_id, a1type=a1type, prop_u_stat=u_stat,
                            prop_p_value=p_value, higher_prop=higher))

    stage = pd.DataFrame.from_records(records)
    if stage.empty:
        return stage
    stage["prop_adj_pvalue"] = significance.benjamini_hochberg(
            stage["prop_p_value"].to_numpy())
    return stage


def site_counts(frame, pairs, case, reference):
    """
    One row per (gene, call type, position) of a retained pair: how many
    individuals of each condition carry it.
    """
    if not pairs:
        return pd.DataFrame(columns=["ref_id", "a1type", "chrom", "pos",
                                     "n_case", "n_ref"])
    wanted = frame.set_index(["ref_id", "a1type"]).index.isin(pairs)
    carriers = (frame.loc[wanted]
                .groupby(["ref_id", "a1type", "chrom", "pos", "condition"])
                ["individualID"].nunique().rename("n").reset_index())
    wide = carriers.pivot_table(
            index=["ref_id", "a1type", "chrom", "pos"],
            columns="condition", values="n", fill_value=0).reset_index()
    for label in (case, reference):
        if label not in wide.columns:
            wide[label] = 0
    wide = wide.rename(columns={case: "n_case", reference: "n_ref"})
    return wide[["ref_id", "a1type", "chrom", "pos", "n_case", "n_ref"]]


#----------------------------------------------------
# writing
#----------------------------------------------------

def check_reloadable(con, overwrite):
    counts = {table: con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
              for table in LOCAL_TABLES}
    if any(counts.values()) and not overwrite:
        raise SystemExit(
                f"tarva analysispipe local-sites: "
                f"{', '.join(f'{t} holds {n:,} row(s)' for t, n in counts.items())}. "
                f"Pass --overwrite to recompute.")


def prepare_tables(con):
    counts = {table: con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
              for table in LOCAL_TABLES}
    if any(counts.values()):
        for table in LOCAL_TABLES:
            con.execute(f"DROP TABLE IF EXISTS {table}")
        con.commit()
        schema.create_schema(con, LOCAL_TABLES)
        schema.report(f"--overwrite: dropped and recreated "
                      f"{', '.join(LOCAL_TABLES)}")


#----------------------------------------------------
# entry point
#----------------------------------------------------

def run(args):
    db_path = Path(args.db)
    if not db_path.is_file():
        raise SystemExit(
                f"tarva analysispipe local-sites: database not found: {db_path}")

    con = schema.connect_for_build(db_path)
    try:
        schema.create_schema(con)
        check_reloadable(con, args.overwrite)

        case, reference, sizes = gene_tests.load_conditions(con, None)
        schema.report(f"conditions: case={case} (n={sizes[case]}), "
                      f"reference={reference} (n={sizes[reference]})")
        schema.report(f"majority threshold: {majority(sizes[case])} of "
                      f"{sizes[case]} {case}, or {majority(sizes[reference])} of "
                      f"{sizes[reference]} {reference}")

        filtered = schema.get_meta(con, "gene_tests.subset_from_global_set")
        if filtered is not None and not int(filtered):
            schema.warn("local-sites: gene_class was built over ALL genes because "
                        "step 1a's global set was empty, so the local stages narrow "
                        "the unfiltered gene set. run_meta records which set was used.")

        frame = load_calls(con, args.min_depth)
        schema.report(f"eligible calls: {len(frame):,} over "
                      f"{frame['ref_id'].nunique():,} gene(s)")

        #---- stage 4
        tests = stage4(frame, case, reference, sizes)
        testable = int(tests["p_value"].notna().sum())
        schema.report(f"stage 4: {len(tests):,} (gene, call type) pair(s), "
                      f"{testable:,} with a defined t-test")

        basis, _ = significance.choose_basis(
                tests, LOCAL_ALPHA, significance.GATING, "stage 4",
                downstream="Stages 5 and 6")
        nominal = tests["p_value"] < NOMINAL_ALPHA
        schema.report(f"stage 4: {int(nominal.sum()):,} pair(s) at nominal p < "
                      f"{NOMINAL_ALPHA}, in {tests.loc[nominal, 'ref_id'].nunique():,} "
                      f"gene(s)")

        if basis is None:
            tests["retained"] = 0
        else:
            tests["retained"] = ((tests[basis] < LOCAL_ALPHA)
                                 & tests["has_majority"]).astype("int64")
        retained = tests.loc[tests["retained"] == 1]
        schema.report(f"stage 4 retained: {len(retained):,} pair(s) in "
                      f"{retained['ref_id'].nunique():,} gene(s) at "
                      f"{basis or 'no qualifying basis'} p < {LOCAL_ALPHA} with a "
                      f"majority of one condition carrying the call type")

        #---- stage 5
        tests["local_dist"] = classify_local(tests, case, reference, sizes)
        for label, count in tests.loc[tests["retained"] == 1,
                                      "local_dist"].value_counts().items():
            schema.report(f"  {label:<14} {count:,}")

        #---- stage 6
        pairs = set(zip(retained["ref_id"], retained["a1type"]))
        stage = stage6(frame, pairs, case, reference)
        tests = tests.merge(stage, on=["ref_id", "a1type"], how="left")
        schema.report(f"stage 6: {int(tests['prop_p_value'].notna().sum()):,} "
                      f"pair(s) with a defined Mann-Whitney")

        #---- final retention
        larger = np.maximum(tests["n_case"], tests["n_ref"])
        smaller = np.minimum(tests["n_case"], tests["n_ref"])
        # The two alternative criteria, kept as columns so the counts below
        # report what actually qualified rather than re-deriving it.
        tests["by_dominance"] = ((larger > 1)
                                 & (larger >= DOMINANCE * smaller)).astype("int64")
        tests["by_proportion"] = (tests["prop_p_value"].notna()
                                  & (tests["prop_p_value"] < NOMINAL_ALPHA)
                                  ).astype("int64")
        tests["final_retained"] = (
                (tests["retained"] == 1)
                & ((tests["by_dominance"] == 1) | (tests["by_proportion"] == 1))
                ).astype("int64")
        final = tests.loc[tests["final_retained"] == 1]
        schema.report(f"final: {len(final):,} pair(s) in "
                      f"{final['ref_id'].nunique():,} gene(s) -- "
                      f"{int(final['by_dominance'].sum()):,} by carrier dominance, "
                      f"{int(final['by_proportion'].sum()):,} by stage-6 p < "
                      f"{NOMINAL_ALPHA} (criteria overlap)")

        sites = site_counts(frame, set(zip(final["ref_id"], final["a1type"])),
                            case, reference)
        schema.report(f"sites: {len(sites):,} altered position(s) in the final "
                      f"pair(s)")

        prepare_tables(con)
        tests.reindex(columns=[
            "ref_id", "a1type", "n_case", "n_ref", "prop_case", "prop_ref",
            "t_stat", "p_value", "adj_pvalue", "retained", "local_dist",
            "prop_u_stat", "prop_p_value", "prop_adj_pvalue", "higher_prop",
            "final_retained"]).to_sql(
                "local_calltype_test", con, if_exists="append", index=False,
                chunksize=gene_tests.WRITE_CHUNK)
        con.commit()
        schema.report(f"local_calltype_test: {len(tests):,} row(s) written")

        sites.to_sql("local_site", con, if_exists="append", index=False,
                     chunksize=gene_tests.WRITE_CHUNK)
        con.commit()
        schema.report(f"local_site: {len(sites):,} row(s) written")

        schema.create_indexes(con, LOCAL_TABLES)
        schema.set_meta(con, "local_sites", {
            "local_sites.case_condition": case,
            "local_sites.reference_condition": reference,
            "local_sites.min_depth": args.min_depth,
            "local_sites.local_alpha": LOCAL_ALPHA,
            "local_sites.nominal_alpha": NOMINAL_ALPHA,
            "local_sites.dominance": DOMINANCE,
            "local_sites.majority_case": majority(sizes[case]),
            "local_sites.majority_reference": majority(sizes[reference]),
            "local_sites.n_pairs": len(tests),
            "local_sites.n_pairs_testable": testable,
            "local_sites.stage4_basis": basis or "none",
            "local_sites.n_genes_nominal": int(
                    tests.loc[tests["p_value"] < NOMINAL_ALPHA, "ref_id"].nunique()),
            "local_sites.n_pairs_retained": len(retained),
            "local_sites.n_pairs_final": len(final),
            "local_sites.n_genes_final": int(final["ref_id"].nunique()),
            "local_sites.n_sites_final": len(sites),
            "local_sites.output_dir": (str(Path(args.output).resolve())
                                       if args.output else ""),
            })
    finally:
        con.close()

    if args.output:
        directory = Path(args.output)
        directory.mkdir(parents=True, exist_ok=True)
        out = tests.copy()
        out["case_condition"] = case
        out["ref_condition"] = reference
        gene_tests.write_csv(out, directory, CALLTYPE_CSV)
        gene_tests.write_csv(sites, directory, SITES_CSV)

    schema.report("")
    schema.report("local-sites complete")
    schema.report(f"  (gene, call type) pairs   {len(tests):,}")
    schema.report(f"  with a defined t-test     {testable:,}")
    schema.report(f"  genes at nominal p < {NOMINAL_ALPHA}  "
                  f"{int(tests.loc[tests['p_value'] < NOMINAL_ALPHA, 'ref_id'].nunique()):,}")
    schema.report(f"  stage 4 retained          {len(retained):,}")
    schema.report(f"  final retained            {len(final):,} pair(s) in "
                  f"{final['ref_id'].nunique():,} gene(s)")
    schema.report(f"  altered positions         {len(sites):,}")
    schema.report(f"Wrote {db_path}")
    return 0
