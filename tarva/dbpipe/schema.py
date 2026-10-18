"""
TARVa database schema

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
    # stored. Every analysis reads this column; none recomputes it. That single
    # change removes the four-way drift recorded in section 4 of the plan.
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

    # complete_in names the condition in which EVERY sample carried at least
    # one eligible call for this (gene, a1type). That completeness is the
    # entry criterion: a pair qualifies only if at least one variant call 
    # is found for all individuals in at least one conditions.
    # Storing it is what makes n_*_observed interpretable
    # without inferring which side was the complete one.
    "gene_test": """
        CREATE TABLE IF NOT EXISTS gene_test (
            ref_id           text NOT NULL,
            a1type           text NOT NULL,
            complete_in      text,
            u_stat           real,
            p_value          real,
            adj_pvalue       real,
            rbc              real,
            n_case_observed  integer,
            n_case_padded    integer,
            n_ref_observed   integer,
            n_ref_padded     integer,
            PRIMARY KEY (ref_id, a1type),
            FOREIGN KEY (ref_id) REFERENCES gene(ref_id)
        ) WITHOUT ROWID;
    """,

    # Separate from `gene` because these are written by `stats`, not by
    # `annotate`, and one writer per table is the rule. Keyed on ref_id alone
    # because edit_dist is a per-gene property -- storing it on gene_test would
    # duplicate it across every call type.
    #
    # Each condition's share of the gene's total, thresholded at 0.75. 
    # The component counts sit beside
    # the verdict for the same reason n_*_observed sits beside a p-value -- so
    # the assignment is checkable rather than asserted.
    "gene_class": """
        CREATE TABLE IF NOT EXISTS gene_class (
            ref_id        text NOT NULL PRIMARY KEY,
            edit_dist     text,
            gene_group    text,
            n_calls_case  integer,
            n_calls_ref   integer,
            sum_ad_case   integer,
            sum_ad_ref    integer,
            FOREIGN KEY (ref_id) REFERENCES gene(ref_id)
        );
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

    # Both keyed on ref_id as PRIMARY KEY, which serves every lookup either
    # one gets; gene_test's (ref_id, a1type) prefix covers per-gene access.
    "gene_test": (),
    "gene_class": (),

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
