"""
analysispipe cascade: count what survived each filter, and the background
rate of each variant call-type.

Every number is a COUNT over `call`, `site` and `gene` before any gene-level
filtering -- no test and no database table. THE STAGES ARE NESTED rather than
independent: a call at a site with no gene always has a NULL call type, since
the type needs the gene's strand, so `classified` is a subset of `in_gene` and
subtracting the two independently would count the same calls twice.

    `tarva analysispipe cascade --db <TARVa.db> [--min-depth 20] [--output <dir>]`

Input:  a database with annotate, load-samples and load-calls done. The two
        gene-level rows are read from run_meta rather than recomputed, since
        `gene-tests` is their single writer; they report as unavailable if it
        has not run. `at_depth` is cross-checked against
        gene_tests.n_eligible_calls, so a cascade run at a different
        --min-depth from the analysis says so instead of being plotted.
Output: filter_cascade.csv (Figure 1) and calltype_rates.csv when --output is
        given, and cascade.* keys in run_meta.
"""

from pathlib import Path

import pandas as pd

from tarva.analysispipe import gene_tests
from tarva.dbpipe import schema

DEFAULT_MIN_DEPTH = 20

CASCADE_CSV = "filter_cascade.csv"
RATES_CSV = "calltype_rates.csv"


def add_arguments(parser):
    parser.add_argument("--db", required=True,
                        help="Database with annotate, load-samples and load-calls done")
    parser.add_argument("--min-depth", type=int, default=DEFAULT_MIN_DEPTH,
                        help=f"The depth stage of the cascade, and the threshold the "
                             f"call-type rates are counted at (default: "
                             f"{DEFAULT_MIN_DEPTH}). Must match the --min-depth "
                             f"`gene-tests` used, which is checked against run_meta")
    parser.add_argument("--output", default=None, metavar="DIR",
                        help="Directory for filter_cascade.csv and calltype_rates.csv. "
                             "Created if absent. Without it the counts are reported to "
                             "stdout and stored in run_meta only")


#----------------------------------------------------
# the cascade
#----------------------------------------------------

def cascade_stages(min_depth):
    """
    (name, description, where, params) per stage, each predicate being the
    previous one plus one term.

    Written as accumulating strings rather than four independent WHERE clauses
    so that the nesting is visible in the source: if a stage's predicate is
    not its predecessor plus a term, the counts are not a cascade.
    """
    in_gene = "s.ref_id IS NOT NULL"
    classified = f"{in_gene} AND c.a1type IS NOT NULL"
    at_depth = f"{classified} AND c.dp >= ?"
    return (
        ("loaded", "every call load-calls wrote", "1 = 1", ()),
        ("in_gene", "at a site resolved to a gene", in_gene, ()),
        ("classified", "with a variant call-type assigned", classified, ()),
        ("at_depth", f"at dp >= {min_depth}", at_depth, (int(min_depth),)),
        )


def count_stage(con, where, params):
    """
    (n_calls, n_sites, n_genes) for one stage.

    One query per stage rather than conditional aggregation in a single pass.
    The table is small enough that four scans cost nothing, and four short
    queries can be read side by side to check the nesting.
    """
    query = f"""
        SELECT COUNT(*),
               COUNT(DISTINCT s.chrom || '_' || s.pos),
               COUNT(DISTINCT s.ref_id)
        FROM call AS c
        JOIN site AS s ON s.chrom = c.chrom AND s.pos = c.pos
        WHERE {where}
    """
    cursor = con.execute(query, params)
    n_calls, n_sites, n_genes = cursor.fetchone()
    return int(n_calls), int(n_sites), int(n_genes)


def build_cascade(con, min_depth):
    """
    One row per stage, with the drop from the stage above it.

    n_dropped is reported rather than left to the reader because a cascade
    figure is drawn from the drops, and computing them from adjacent rows is
    where an off-by-one stage ordering would hide.

    A NOTE ON n_genes. COUNT(DISTINCT s.ref_id) ignores NULLs, so the `loaded`
    stage's gene count already excludes geneless sites and necessarily equals
    `in_gene`'s. That is not the query failing to filter: the calls dropped
    between those stages belong to no gene, so dropping them cannot remove
    one. Equal gene counts across the first three stages is the correct and
    informative answer -- only the depth stage costs genes.
    """
    schema.require_table(con, "call", "tarva dbpipe load-calls")
    schema.require_table(con, "gene", "tarva dbpipe annotate")

    records, previous = [], None
    for name, description, where, params in cascade_stages(min_depth):
        n_calls, n_sites, n_genes = count_stage(con, where, params)
        records.append({
            "stage": name, "description": description,
            "n_calls": n_calls, "n_sites": n_sites, "n_genes": n_genes,
            "n_calls_dropped": 0 if previous is None else previous - n_calls,
            })
        previous = n_calls
    return pd.DataFrame.from_records(records)


#----------------------------------------------------
# background call-type rates
#----------------------------------------------------

def calltype_rates(con, min_depth):
    """
    One row per (call type, condition): calls, distinct sites, individuals.

    Both counts are reported: the total variant calls across all samples and
    all genes, and the distinct positions per call type. Ordered by site count
    descending, so the rows read as a frequency ordering.

    n_calls IS additive across rows and sums to the cascade's at_depth stage.
    n_sites, n_individuals and n_genes are NOT: each is a distinct count
    within its own (call type, condition) group, so a position called in both
    conditions is counted in both rows, and a position where two individuals
    carry different alternate alleles is counted under each call type. Summing
    those columns overcounts. Use n_calls for a total and read the others per
    row.
    """
    table = schema.INDIVIDUAL_TABLE
    schema.require_table(con, table, "tarva dbpipe load-samples")
    id_column, condition_column = schema.individual_columns(con)
    query = f"""
        SELECT c.a1type AS a1type,
               p."{condition_column}" AS condition,
               COUNT(*) AS n_calls,
               COUNT(DISTINCT s.chrom || '_' || s.pos) AS n_sites,
               COUNT(DISTINCT c.individualID) AS n_individuals,
               COUNT(DISTINCT s.ref_id) AS n_genes
        FROM call AS c
        JOIN site AS s ON s.chrom = c.chrom AND s.pos = c.pos
        JOIN {table} AS p ON p."{id_column}" = c.individualID
        WHERE s.ref_id IS NOT NULL
          AND c.a1type IS NOT NULL
          AND c.dp >= ?
        GROUP BY c.a1type, p."{condition_column}"
        ORDER BY n_sites DESC, c.a1type
    """
    return pd.read_sql_query(query, con, params=(int(min_depth),))


#----------------------------------------------------
# what gene-tests already recorded
#----------------------------------------------------

DOWNSTREAM_KEYS = (
    ("gene_tests.n_genes_considered", "genes with an eligible call"),
    ("gene_tests.n_genes_complete", "genes passing the entry criterion"),
    ("gene_tests.n_genes_global_set", "genes in the global gene set"),
    ("gene_tests.n_calls_in_subsets", "calls in the genes steps 2a/2b used"),
    )


def downstream(con):
    """
    [(label, value)] for the gene-level stages, read from run_meta.

    Not recomputed: `gene-tests` owns the entry criterion and the global gene
    set, and a second implementation of completeness here would be free to
    disagree with the one the analysis actually used.
    """
    rows = []
    for key, label in DOWNSTREAM_KEYS:
        rows.append((label, schema.get_meta(con, key)))
    return rows


def check_depth(con, min_depth, at_depth_calls):
    """
    Warn when this command's --min-depth disagrees with the one gene-tests used.

    The cascade only describes the analysis's own filtering if the two match.
    """
    recorded = schema.get_meta(con, "gene_tests.min_depth")
    if recorded is not None and int(recorded) != int(min_depth):
        schema.warn(f"cascade: --min-depth {min_depth} does not match the "
                    f"{recorded} gene-tests used. The depth stage below describes a "
                    f"different filtering than the analysis performed; re-run with "
                    f"--min-depth {recorded} for a cascade that matches it.")
        return
    eligible = schema.get_meta(con, "gene_tests.n_eligible_calls")
    if eligible is not None and int(eligible) != at_depth_calls:
        schema.warn(f"cascade: the depth stage holds {at_depth_calls:,} call(s) but "
                    f"gene-tests recorded {int(eligible):,} eligible. The two should be "
                    f"identical at the same --min-depth; one of the tables has changed "
                    f"since that run.")


#----------------------------------------------------
# entry point
#----------------------------------------------------

def run(args):
    db_path = Path(args.db)
    if not db_path.is_file():
        raise SystemExit(f"tarva analysispipe cascade: database not found: {db_path}")

    con = schema.connect_for_build(db_path)
    try:
        stages = build_cascade(con, args.min_depth)
        rates = calltype_rates(con, args.min_depth)

        at_depth_calls = int(stages.loc[stages["stage"] == "at_depth",
                                        "n_calls"].iloc[0])
        check_depth(con, args.min_depth, at_depth_calls)

        schema.report(f"filtering cascade (dp >= {args.min_depth}):")
        for row in stages.itertuples():
            dropped = ("" if not row.n_calls_dropped
                       else f"  -{row.n_calls_dropped:,}")
            schema.report(f"  {row.stage:<11} {row.n_calls:>10,} call(s)  "
                          f"{row.n_sites:>9,} site(s)  {row.n_genes:>7,} gene(s)"
                          f"{dropped}")
            schema.report(f"              {row.description}")

        schema.report("")
        schema.report("gene-level stages, from run_meta:")
        for label, value in downstream(con):
            shown = "not recorded -- run gene-tests" if value is None else f"{int(value):,}"
            schema.report(f"  {label:<40} {shown}")

        schema.report("")
        schema.report(f"call-type rates over all genes (dp >= {args.min_depth}):")
        for a1type, block in rates.groupby("a1type", sort=False):
            parts = "  ".join(f"{row.condition} {row.n_calls:,}/{row.n_sites:,}"
                              for row in block.itertuples())
            schema.report(f"  {a1type:<8} {parts}")
        schema.report("  (condition calls/sites)")

        schema.set_meta(con, "cascade", {
            "cascade.min_depth": args.min_depth,
            "cascade.n_calls_loaded": int(stages["n_calls"].iloc[0]),
            "cascade.n_calls_in_gene": int(
                    stages.loc[stages["stage"] == "in_gene", "n_calls"].iloc[0]),
            "cascade.n_genes_in_gene": int(
                    stages.loc[stages["stage"] == "in_gene", "n_genes"].iloc[0]),
            "cascade.n_calls_classified": int(
                    stages.loc[stages["stage"] == "classified", "n_calls"].iloc[0]),
            "cascade.n_calls_at_depth": at_depth_calls,
            "cascade.n_calltypes": int(rates["a1type"].nunique()),
            "cascade.output_dir": (str(Path(args.output).resolve())
                                   if args.output else ""),
            })
    finally:
        con.close()

    if args.output:
        directory = Path(args.output)
        directory.mkdir(parents=True, exist_ok=True)
        gene_tests.write_csv(stages, directory, CASCADE_CSV)
        gene_tests.write_csv(rates, directory, RATES_CSV)

    schema.report("")
    schema.report("cascade complete")
    for row in stages.itertuples():
        schema.report(f"  {row.stage:<11} {row.n_calls:>10,} call(s)  "
                      f"{row.n_genes:>7,} gene(s)")
    schema.report(f"  call types                {rates['a1type'].nunique()}")
    schema.report(f"Wrote {db_path}")
    return 0
