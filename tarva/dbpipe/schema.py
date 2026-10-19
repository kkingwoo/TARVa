"""
The TARVa database schema: every table, its indexes, and the helpers the
build and analysis steps share.

One step is the single writer of any given table, and each step records its
own thresholds and input paths in run_meta under its own key prefix.

Input:  nothing -- this module is imported, not run.
Output: table and index definitions, connection helpers
        (`connect_for_build`, `connect_read_only`), schema creation
        (`create_schema`, `create_indexes`) and run_meta access
        (`set_meta`, `get_meta`).
"""

import sqlite3 as sql
import sys

from datetime import datetime
from pathlib import Path

SCHEMA_VERSION = "1"

#----------------------------------------------------
# Create tables
#----------------------------------------------------

TABLES = {

    # NOTE:`individual` is NOT declared here. It is the one table whose
    # columns are a property of the DATASET rather than of this design: it is
    # the metadata CSV itself, written verbatim by `load-samples` with pandas
    # to_sql, so every column name transfers and nothing is lost by omission.

    # ref_id is the UNVERSIONED reference gene identifier with any GFF3 ID
    # namespace removed, or the assembly locus id (MSTRG.203) for a locus with
    # no reference gene at all.
    "gene": """
        CREATE TABLE IF NOT EXISTS gene (
            ref_id          text    NOT NULL PRIMARY KEY,
            ref_id_version  text,
            gene_name       text,
            biotype         text,
            novel           integer NOT NULL,
            chrom           text    NOT NULL,
            strand          text    NOT NULL,
            gene_start      integer NOT NULL,
            gene_end        integer NOT NULL,
            gene_span       integer NOT NULL,
            n_transcripts   integer NOT NULL
        );
    """,

    # ref_id is NULL when the transcript belongs to a chimeric locus and has no
    # unambiguous parent gene (a novel transcript in a locus carrying two or
    # more reference genes). 
    "transcript": """
        CREATE TABLE IF NOT EXISTS transcript (
            tid              text    NOT NULL PRIMARY KEY,
            tid_version      text,
            ref_id           text,
            locus_id         text    NOT NULL,
            chrom            text    NOT NULL,
            strand           text    NOT NULL,
            trans_start      integer NOT NULL,
            trans_end        integer NOT NULL,
            trans_len        integer NOT NULL,
            transcript_type  text,
            tsl              integer,
            ccds             text,
            FOREIGN KEY (ref_id) REFERENCES gene(ref_id)
        );
    """,

    # exon_num is order along the transcript, 1 = 5'-most. trans_offset is the
    # number of transcript bases preceding this exon in that order, so a
    # transcript coordinate is one addition:
    #   +  trans_pos = trans_offset + (pos - exon_start) + 1
    #   -  trans_pos = trans_offset + (exon_end - pos) + 1
    # Coordinates are the GTF's own 1-based inclusive values, so an exon's
    # length is (exon_end - exon_start + 1).
    "exon": """
        CREATE TABLE IF NOT EXISTS exon (
            tid           text    NOT NULL,
            exon_num      integer NOT NULL,
            exon_start    integer NOT NULL,
            exon_end      integer NOT NULL,
            trans_offset  integer NOT NULL,
            PRIMARY KEY (tid, exon_num),
            FOREIGN KEY (tid) REFERENCES transcript(tid)
        ) WITHOUT ROWID;
    """,

    # One row per distinct genomic position with a call. ref_base is
    # the GENOMIC reference nucleotide.
    #
    # ref_base is the single genomic base AT pos -- the first character of the
    # VCF REF, which is a true function of the coordinate. It is NOT the REF
    # allele: an indel's REF spans the affected bases, so two samples at one
    # position legitimately carry REF=TGG and REF=T (a deletion and an
    # insertion sharing an anchor). The full per-call allele lives in
    # call.ref.
    #
    # ref_id NULL means no gene span covers this position, and strand is then
    # NULL too -- strand is the covering GENE's strand, for every variant
    # class alike.
    "site": """
        CREATE TABLE IF NOT EXISTS site (
            chrom      text    NOT NULL,
            pos        integer NOT NULL,
            strand     text,
            ref_base   text    NOT NULL,
            ref_id     text,
            PRIMARY KEY (chrom, pos),
            FOREIGN KEY (ref_id) REFERENCES gene(ref_id)
        ) WITHOUT ROWID;
    """,

    # ref and alt1 are this call's own VCF alleles.
    #
    # a1_ad, dp, ref and alt1 are NOT NULL: they come straight from the VCF,
    # and a variant call inherently carries all four.
    #
    # a1type is computed once, at load, from (ref, alt1, site.strand) and
    # stored. Every analysis reads this column and none recomputes it, so the
    # call type cannot differ between steps.
    "call": """
        CREATE TABLE IF NOT EXISTS call (
            chrom   text    NOT NULL,
            pos     integer NOT NULL,
            individualID  text NOT NULL,
            ref     text    NOT NULL,
            alt1    text    NOT NULL,
            a1_ad   integer NOT NULL,
            dp      integer NOT NULL,
            a1type  text,
            filter  text    NOT NULL,
            PRIMARY KEY (chrom, pos, individualID),
            FOREIGN KEY (chrom, pos) REFERENCES site(chrom, pos),
            FOREIGN KEY (individualID) REFERENCES individual(individualID)
        ) WITHOUT ROWID;
    """,

    "expression": """
        CREATE TABLE IF NOT EXISTS expression (
            individualID  text NOT NULL,
            tid   text NOT NULL,
            FPKM  real,
            TPM   real,
            PRIMARY KEY (individualID, tid),
            FOREIGN KEY (individualID) REFERENCES individual(individualID),
            FOREIGN KEY (tid) REFERENCES transcript(tid)
        ) WITHOUT ROWID;
    """,

    # The per-gene test: ONE ROW PER GENE, every call type pooled. Each individual
    # contributes a single number -- all of their unweighted ad/dp proportions
    # in that gene, across every position AND every call type, summed -- and
    # the arms are those per-individual sums, zero-padded to cohort size.
    #
    # Per-call-type quantities live in gene_calltype, not here: the unit of
    # this test is the gene, and storing a call type against a per-gene
    # p-value would invite reading the two as if they belonged together.
    #
    # complete_in names the condition in which EVERY individual carried at
    # least one eligible call in this gene -- the entry criterion. Storing it
    # is what makes n_*_observed interpretable without inferring which side
    # was the complete one.
    #
    # n_sites_* are raw SITE counts per condition, not individual counts,
    # tracked beside the per-individual sums the test actually uses.
    "gene_test": """
        CREATE TABLE IF NOT EXISTS gene_test (
            ref_id           text NOT NULL PRIMARY KEY,
            complete_in      text,
            u_stat           real,
            p_value          real,
            adj_pvalue       real,
            rbc              real,
            n_case_observed  integer,
            n_case_padded    integer,
            n_ref_observed   integer,
            n_ref_padded     integer,
            n_sites_case     integer,
            n_sites_ref      integer,
            FOREIGN KEY (ref_id) REFERENCES gene(ref_id)
        );
    """,

    # Per (gene, call type) aggregates. NO TEST HERE -- local stage 4 runs its
    # own t-test and is the only consumer. This table exists because that test
    # needs quantities at this unit and nothing else in the pipeline does.
    #
    # sum_ad_* / sum_dp_* give the aggregated proportion sum_ad / sum_dp per
    # gene per call type per condition: allele depths and supporting read
    # depths summed over every alternate-allele site and every individual of
    # that condition.
    #
    # n_indiv_* count DISTINCT INDIVIDUALS carrying that call type in that
    # gene, which is what stage 4's retention rule tests -- at least half the
    # individuals of one condition. n_sites_* count rows. The two differ
    # whenever one individual carries several sites, which is the usual case,
    # so both are stored rather than one inferred from the other.
    "gene_calltype": """
        CREATE TABLE IF NOT EXISTS gene_calltype (
            ref_id        text NOT NULL,
            a1type        text NOT NULL,
            n_sites_case  integer,
            n_sites_ref   integer,
            n_indiv_case  integer,
            n_indiv_ref   integer,
            sum_ad_case   integer,
            sum_dp_case   integer,
            sum_ad_ref    integer,
            sum_dp_ref    integer,
            PRIMARY KEY (ref_id, a1type),
            FOREIGN KEY (ref_id) REFERENCES gene(ref_id)
        ) WITHOUT ROWID;
    """,

    # Separate from `gene` because these are written by `gene-tests`, not by
    # `annotate`, and one writer per table is the rule.
    #
    # TWO three-category classifications of the same gene, at two thresholds:
    #
    #   gene_group  the condition-enriched subset assignment. A condition
    #               holding >= 75% of the gene's alternate-allele COUNTS gives
    #               that condition's *_enriched subset; anything between 25%
    #               and 75% gives Shared.
    #   edit_dist   the per-gene distribution label the site-uniqueness step
    #               groups by. Whichever condition holds MORE THAN HALF gives
    #               that condition's *_enriched label; an exact tie gives Both.
    #
    # Both name two of their groups *_enriched; only the threshold and the
    # third label differ (Shared at 75%, Both at a tie).
    #
    # The component counts sit beside the verdict for the same reason
    # n_*_observed sits beside a p-value.
    "gene_class": """
        CREATE TABLE IF NOT EXISTS gene_class (
            ref_id        text NOT NULL PRIMARY KEY,
            edit_dist     text,
            gene_group    text,
            n_calls_case  integer,
            n_calls_ref   integer,
            FOREIGN KEY (ref_id) REFERENCES gene(ref_id)
        );
    """,

    # One row per gene: Spearman between per-gene TPM and variant count across
    # the individuals of each condition, then Fisher z on the difference
    # between the two correlations.
    #
    # case/ref rather than the dataset's own labels, as everywhere else, with
    # the labels in run_meta.
    #
    # NULL rho and p_value mean the condition failed the entry rule for this
    # test (fewer than MIN_INDIVIDUALS individuals, or no variance in TPM or
    # in variant count), which is a different statement from a correlation of
    # zero. 
    #
    # Both diff_pvalue and diff_adj_pvalue are stored, since the figure
    # outlines points on the nominal value while the adjusted one decides
    # significance.
    "expression_test": """
        CREATE TABLE IF NOT EXISTS expression_test (
            ref_id           text NOT NULL PRIMARY KEY,
            gene_group       text,
            novel            integer,
            n_case           integer,
            case_rho         real,
            case_pvalue      real,
            n_ref            integer,
            ref_rho          real,
            ref_pvalue       real,
            diff_pvalue      real,
            diff_adj_pvalue  real,
            FOREIGN KEY (ref_id) REFERENCES gene(ref_id)
        );
    """,

    # One row per (edit_dist subset, call type). A 2x2 of condition against
    # site uniqueness, tested with Fisher's exact.
    #
    #                 unique            common
    #   case     n_case_unique     n_case_common
    #   ref      n_ref_unique      n_ref_common
    #
    # odds_ratio is (case_unique * ref_common) / (case_common * ref_unique),
    # so ABOVE 1 means unique sites are relatively more likely in the case
    # condition. The four counts are distinct POSITIONS, not calls.
    #
    # p_value comes from Fisher's exact on the RAW counts -- an exact test
    # needs no zero correction. odds_ratio and the interval come from counts
    # with Haldane-Anscombe applied, because a zero cell makes the ratio
    # infinite or undefined. Both are stored so neither has to be recomputed
    # from the other's inputs.
    "uniqueness_test": """
        CREATE TABLE IF NOT EXISTS uniqueness_test (
            edit_dist      text NOT NULL,
            a1type         text NOT NULL,
            n_case_unique  integer,
            n_case_common  integer,
            n_ref_unique   integer,
            n_ref_common   integer,
            odds_ratio     real,
            ci_low         real,
            ci_high        real,
            p_value        real,
            adj_pvalue     real,
            basis          text,
            PRIMARY KEY (edit_dist, a1type)
        ) WITHOUT ROWID;
    """,

    # One test per call type, across genes. gene_set names which set of genes
    # the test ran over. One observation is one INDIVIDUAL -- their weighted
    # proportions of that call type summed over every position and every gene
    # in the set -- so the arms are zero-padded to cohort size as gene_test's
    # are, and the observed/padded split is recorded for the same reason: a
    # result resting mostly on imputed zeros should be visible without a
    # re-run.
    #
    # `basis` names the p-value column this step's significance was read from
    # -- adj_pvalue, p_value, or empty when neither reached alpha in a
    # majority of tests. Stored because the call-type figure must mark the
    # column the step actually used, which it cannot infer from the values.
    "calltype_test": """
        CREATE TABLE IF NOT EXISTS calltype_test (
            gene_set         text NOT NULL,
            a1type           text NOT NULL,
            u_stat           real,
            p_value          real,
            adj_pvalue       real,
            basis            text,
            rbc              real,
            n_case_observed  integer,
            n_case_padded    integer,
            n_ref_observed   integer,
            n_ref_padded     integer,
            PRIMARY KEY (gene_set, a1type)
        ) WITHOUT ROWID;
    """,

    # The local level, stages 4-6: one row per (gene, call type).
    #
    # stage 4  prop_* is that condition's aggregated proportion, sum(ad) over
    #          sum(dp) across every alternate-allele site and every CARRIER of
    #          the condition. 
    # stage 5  n_case / n_ref count DISTINCT INDIVIDUALS carrying that call
    #          type in that gene, and local_dist classifies the pair from them.
    # stage 6  prop_p_value is a Mann-Whitney over the length-weighted
    #          per-site proportions, pooled across the gene's positions and
    #          split by condition, with higher_prop naming the condition whose
    #          mean was larger.
    #
    # retained is stage 4's: p below its alpha AND at least half the
    # individuals of one condition carrying the call type. final_retained adds
    # stage 6's rule -- one condition's carrier count above 1 and at least
    # half again the other's, or a nominal stage-6 p below its alpha.
    "local_calltype_test": """
        CREATE TABLE IF NOT EXISTS local_calltype_test (
            ref_id           text NOT NULL,
            a1type           text NOT NULL,
            n_case           integer,
            n_ref            integer,
            prop_case        real,
            prop_ref         real,
            t_stat           real,
            p_value          real,
            adj_pvalue       real,
            retained         integer,
            local_dist       text,
            prop_u_stat      real,
            prop_p_value     real,
            prop_adj_pvalue  real,
            higher_prop      text,
            final_retained   integer,
            PRIMARY KEY (ref_id, a1type),
            FOREIGN KEY (ref_id) REFERENCES gene(ref_id)
        ) WITHOUT ROWID;
    """,

    # One row per altered position of a retained (gene, call type), with how
    # many individuals of each condition carry it.
    #
    # Genomic chrom/pos, as every other table keys positions. A
    # transcript-relative coordinate is computable from `exon.trans_offset`
    # if it is wanted.
    "local_site": """
        CREATE TABLE IF NOT EXISTS local_site (
            ref_id   text    NOT NULL,
            a1type   text    NOT NULL,
            chrom    text    NOT NULL,
            pos      integer NOT NULL,
            n_case   integer,
            n_ref    integer,
            PRIMARY KEY (ref_id, a1type, chrom, pos),
            FOREIGN KEY (chrom, pos) REFERENCES site(chrom, pos)
        ) WITHOUT ROWID;
    """,

    # Thresholds, cohort sizes, tool versions and input paths.
    "run_meta": """
        CREATE TABLE IF NOT EXISTS run_meta (
            key         text NOT NULL PRIMARY KEY,
            value       text,
            written_by  text,
            written_at  text
        );
    """,
}

# Indices assigned once tables are built.
INDEXES = {

    # None here: `individual` carries the metadata CSV's own column names, so
    # the name of its identifier column is not known until load time. That
    # step creates the UNIQUE index itself, against whatever --id-column
    # named. ANALYZE still runs.
    "individual": (),

    "gene": (
        # the competitor lookup -- every gene whose span covers a position
        ("gene_locus_idx",  "CREATE INDEX IF NOT EXISTS gene_locus_idx "
                            "ON gene(chrom, gene_start, gene_end);"),
        ("gene_name_idx",   "CREATE INDEX IF NOT EXISTS gene_name_idx "
                            "ON gene(gene_name);"),
        ),

    "transcript": (
        ("trans_gene_idx",  "CREATE INDEX IF NOT EXISTS trans_gene_idx "
                            "ON transcript(ref_id);"),
        ("trans_locus_idx", "CREATE INDEX IF NOT EXISTS trans_locus_idx "
                            "ON transcript(locus_id);"),
        ("trans_span_idx",  "CREATE INDEX IF NOT EXISTS trans_span_idx "
                            "ON transcript(chrom, trans_start, trans_end);"),
        ),

    # None: PRIMARY KEY (tid, exon_num) on a WITHOUT ROWID table IS the index,
    # and every lookup is a tid prefix scan of it. ANALYZE still runs.
    "exon": (),

    # per-gene site selection, the entry point of nearly every analysis query
    "site": (
        ("site_gene_idx",   "CREATE INDEX IF NOT EXISTS site_gene_idx "
                            "ON site(ref_id, pos);"),
        ),

    # per-sample access; the PK already serves (chrom, pos) prefix scans
    "call": (
        ("call_individual_idx",
                            "CREATE INDEX IF NOT EXISTS call_individual_idx "
                            "ON call(individualID);"),
        ("call_type_idx",   "CREATE INDEX IF NOT EXISTS call_type_idx ON call(a1type);"),
        ),

    "expression": (
        ("expr_tid_idx",    "CREATE INDEX IF NOT EXISTS expr_tid_idx "
                            "ON expression(tid);"),
        ),

    # All three keyed on ref_id as PRIMARY KEY or as its leading column, which
    # serves every lookup they get; gene_calltype's (ref_id, a1type) prefix
    # covers per-gene access as well as the exact key.
    "gene_test": (),
    "gene_class": (),
    "gene_calltype": (),

    # Keyed on ref_id; the only other access is a scan for the figure.
    "expression_test": (),

    # None: PRIMARY KEY (edit_dist, a1type) on a WITHOUT ROWID table is the
    # index, and the table holds three subsets times fourteen call types.
    "uniqueness_test": (),

    # (ref_id, a1type) leads both keys, which serves per-gene access as well
    # as the exact key; local_site also gets a position index, since Table 3
    # is read by position as often as by gene.
    "local_calltype_test": (),
    "local_site": (
        ("local_site_pos_idx",
                            "CREATE INDEX IF NOT EXISTS local_site_pos_idx "
                            "ON local_site(chrom, pos);"),
        ),

    # None: PRIMARY KEY (gene_set, a1type) on a WITHOUT ROWID table is the
    # index, and the table holds one row per call type per gene set -- tens of
    # rows, not millions. ANALYZE still runs.
    "calltype_test": (),

    "run_meta": (),
    }


#----------------------------------------------------
# connections
#----------------------------------------------------

def connect_for_build(db_path):
    """
    A writer connection tuned for bulk loading.

    """
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = sql.connect(str(db_path))
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA temp_store=MEMORY")
    con.execute("PRAGMA cache_size=-262144")   
    con.execute("PRAGMA foreign_keys=OFF")  
    return con


def connect_read_only(db_path, busy_timeout_ms=30000):
    """
    A reader connection. One per thread -- never shared across threads.

    """
    db_path = Path(db_path)
    if not db_path.is_file():
        raise SystemExit(f"tarva dbpipe: database not found: {db_path}")
    con = sql.connect(f"{db_path.resolve().as_uri()}?mode=ro", uri=True)
    con.execute(f"PRAGMA busy_timeout={int(busy_timeout_ms)}")
    return con


def create_schema(con, tables=None):
    """Create the declared tables. Idempotent."""
    wanted = tables or TABLES.keys()
    for name in wanted:
        con.execute(TABLES[name])
    con.commit()


def create_indexes(con, tables=None):
    """
    Build the indexes for `tables` and ANALYZE them. Called by each step once
    its own rows are inserted.
    """
    wanted = list(tables) if tables else list(INDEXES)
    unknown = [table for table in wanted if table not in INDEXES]
    if unknown:
        raise KeyError(f"dbpipe schema: no index set declared for table(s) "
                       f"{', '.join(unknown)} -- add them to INDEXES")

    built = []
    for table in wanted:
        for name, statement in INDEXES[table]:
            con.execute(statement)
            built.append(name)
        con.execute(f"ANALYZE {table}")
    con.commit()
    report(f"indexed {', '.join(wanted)}"
           + (f" -- {', '.join(built)}" if built else " -- primary keys only"))
    return built


#----------------------------------------------------
# run_meta
#----------------------------------------------------

def set_meta(con, written_by, values):
    """
    Record provenance. `values` is a {key: value} mapping; keys are namespaced
    by step (e.g. "annotate.gtf") so each step owns its own key space and the
    one-writer rule holds for run_meta too.
    """
    stamp = datetime.now().isoformat(timespec="seconds")
    rows = [(str(k), None if v is None else str(v), written_by, stamp)
            for k, v in values.items()]
    fill_in = ("INSERT INTO run_meta (key, value, written_by, written_at) "
               "VALUES (?,?,?,?) "
               "ON CONFLICT(key) DO UPDATE SET value=excluded.value, "
               "written_by=excluded.written_by, written_at=excluded.written_at")
    con.executemany(fill_in, rows)
    con.commit()


def get_meta(con, key, default=None):
    row = con.execute("SELECT value FROM run_meta WHERE key = ?", (key,)).fetchone()
    return default if row is None else row[0]


INDIVIDUAL_TABLE = "individual"


def individual_columns(con):
    """
    (id_column, condition_column) for the `individual` table.

    """
    id_column = get_meta(con, "load_samples.id_column")
    condition_column = get_meta(con, "load_samples.condition_column")
    if not id_column or not condition_column:
        raise SystemExit(
                "tarva dbpipe: run_meta has no load_samples.id_column or "
                "load_samples.condition_column, so which column of `individual` holds "
                "the identifier is unknown -- run `dbpipe load-samples` first.")
    return id_column, condition_column


def tarva_version():
    try:
        import tarva
        return getattr(tarva, "__version__", "unknown")
    except Exception:
        return "unknown"


def require_table(con, table, step):
    """Fail with a useful message when a prerequisite step has not run."""
    found = con.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?", (table,)
            ).fetchone()
    if found is None:
        raise SystemExit(f"tarva dbpipe: table {table} does not exist -- run `{step}` first")
    count = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    if not count:
        raise SystemExit(f"tarva dbpipe: table {table} is empty -- run `{step}` first")
    return count


def report(message):
    print(message, flush=True)


def warn(message):
    print(f"WARNING: {message}", file=sys.stderr, flush=True)
