"""
analysispipe external-checks: are TARVa's A-G calls known RNA editing sites?

    `tarva analysispipe external-checks --db <TARVa.db>
        --redi <REDIportal.txt> --editome <ADeditome.txt> --output <results_dir>
        [--external-db PATH] [--chunk-size 500000] [--reload]
        [--editome-columns LIST]
        [--table All_Info_TopGenes_tab] [--a1type A-I]
        [--chrom-column chrom] [--pos-column poss]
        [--hpc] [--conda-env tarva] [--directives-dict "{...}"]`

    An A-G call is the signature of ADAR deamination, but nothing in the TARVa
    build distinguishes a genuine editing site from a mis-called variant. Two
    external catalogues answer different halves of that: REDIportal says the
    position is edited in normal tissue, the AD editome says it was reported as
    edited in Alzheimer's. They are checked independently and reported in
    separate columns, since a position can be in both, one, or neither, and
    which of those it is carries the interpretation.

    --a1type selects the call-type label to check, and it must match what the
    database actually holds for an A->G substitution. The label is already
    resolved for strand, so a minus-strand T>C is included; indel labels are
    not.

    The two source files are tens of millions of rows, so they are streamed in
    chunks into their own database (ExternalChecks.db, WAL). Only the distinct
    Region/Position pairs of the normal-editing catalogue are kept -- presence
    is all the check asks of it -- while the disease catalogue keeps its full
    rows, since ADeditom_ID is reported. TARVa.db is never written to: it is
    ATTACHed read-only and the join runs inside SQLite, so only the result --
    one row per distinct position -- is ever materialised in pandas.

    Loading is skipped when a table is already populated, so re-running the
    check against a different --db costs nothing. Pass --reload to rebuild.

    Matching is POSITIONAL: chrom/pos only, no base or strand comparison. A
    REDIportal row on the minus strand is recorded as T->C in genome
    coordinates while being A->G on the transcript, so requiring the bases to
    agree would silently drop half the catalogue.

    Input:  the TARVa database, plus the two catalogue files as --redi and
            --editome. Loading is skipped when a table is already populated,
            so re-running against a different --db costs nothing; --reload
            rebuilds.
    Output: into --output -- ExternalChecks.db (redi_tab, editome_tab,
            ag_external_check_tab), 09_ag_external_check.csv,
            09_external_check_summary.json
"""

import json
import sqlite3 as sql
import sys

from pathlib import Path

import pandas as pd

from tarva.analysispipe import db as tarva_db
from tarva.slurm_utils import parse_directives_dict, split_directives, submit_job, write_sbatch_script

DEFAULT_EXTERNAL_DB = "ExternalChecks.db"
DEFAULT_CHUNK_SIZE = 500000
DEFAULT_A1TYPE = "A-I"
DEFAULT_CHROM_COLUMN = "chrom"
DEFAULT_POS_COLUMN = "poss"

REDI_TABLE = "redi_tab"
REDI_CHROM, REDI_POS = "Region", "Position"

EDITOME_TABLE = "editome_tab"
EDITOME_CHROM, EDITOME_POS, EDITOME_ID = "chromosome", "editing_position", "ADeditom_ID"

RESULT_TABLE = "ag_external_check_tab"


def _csv_list(raw):
    return [v.strip() for v in str(raw).split(",") if v.strip()] if raw else None


#----------------------------------------------------
# ExternalChecks.db -- loading the two source files
#----------------------------------------------------

def connect_external(path):
    """
    WAL so the load's writes never block a reader, and check_same_thread=False
    because each table is written exactly once, in one call, from one place --
    the connection is then only ever read from.
    """
    con = sql.connect(str(path), check_same_thread=False, uri=True)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA synchronous=NORMAL")
    return con


def _row_count(con, table):
    db = con.cursor()
    check_tab = "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ?"
    if db.execute(check_tab, (table,)).fetchone() is None:
        return None
    return db.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]


def load_redi(con, path, chunk_size, reload_=False):
    """
    Load only the DISTINCT (Region, Position) pairs of the normal-editing
    catalogue.

    The check asks exactly one thing of this file -- is the position in it --
    so none of its other 23 columns are needed, and dropping them takes the
    table from the size of the source file to a fraction of it.

    Deduplication is done by the schema, not by pandas: a chunk can only be
    made unique within itself, whereas (Region, Position) must be unique across
    the whole file. A WITHOUT ROWID primary key with INSERT OR IGNORE gets that
    for free and leaves the table needing no separate index.
    """
    existing = _row_count(con, REDI_TABLE)
    if existing and not reload_:
        print(f"{REDI_TABLE}: {existing:,} unique position(s) already loaded -- skipping "
              f"(pass --reload to rebuild)")
        return existing

    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"tarva analysispipe: --redi file not found: {path}")

    print(f"Loading distinct {REDI_CHROM},{REDI_POS} pairs from {path.name} into "
          f"{REDI_TABLE} (chunks of {chunk_size:,} rows)")

    db = con.cursor()
    make_tab = f"""
        CREATE TABLE {REDI_TABLE} (
            {REDI_CHROM} varchar,
            {REDI_POS} integer,
            PRIMARY KEY ({REDI_CHROM},{REDI_POS})) WITHOUT ROWID;
    """
    db.execute(f"DROP TABLE IF EXISTS {REDI_TABLE}")
    db.execute(make_tab)
    con.commit()

    fill_in = f"INSERT OR IGNORE INTO {REDI_TABLE} ({REDI_CHROM},{REDI_POS}) VALUES(?,?)"
    reader = pd.read_csv(
            path, sep="\t", dtype=str, chunksize=chunk_size, usecols=[REDI_CHROM, REDI_POS])

    read, dropped = 0, 0
    try:
        for n, chunk in enumerate(reader, start=1):
            chunk.columns = [str(c).lstrip("#").strip() for c in chunk.columns]
            read += len(chunk)

            chroms = chunk[REDI_CHROM].astype(str).str.strip()
            positions = pd.to_numeric(chunk[REDI_POS], errors="coerce")
            keep = positions.notna()
            dropped += int((~keep).sum())

            pairs = pd.DataFrame({
                "chrom": chroms[keep].values,
                "pos": positions[keep].astype("int64").values,
                }).drop_duplicates()
            # .tolist() rather than itertuples: sqlite3 cannot bind numpy.int64.
            db.executemany(fill_in, list(zip(pairs["chrom"].tolist(), pairs["pos"].tolist())))
            con.commit()
            print(f"  {REDI_TABLE}: chunk {n} -- {read:,} rows read", flush=True)
    except ValueError as exc:
        raise SystemExit(
                f"tarva analysispipe: could not read {REDI_CHROM} and {REDI_POS} from "
                f"{path}: {exc}"
                )

    total = _row_count(con, REDI_TABLE)
    print(f"{REDI_TABLE}: {total:,} unique position(s) from {read:,} rows")
    if dropped:
        print(f"WARNING: {dropped:,} row(s) of {path.name} had a non-numeric {REDI_POS} "
              f"and were skipped", file=sys.stderr)
    return total


def load_source(con, path, table, chrom_col, pos_col, chunk_size,
                columns=None, reload_=False, label=None):
    """
    Stream one tab-delimited catalogue into `table`, then index it on
    (chrom, pos) -- the only access pattern the check has, and the difference
    between a scan of tens of millions of rows per position and a lookup.

    Positions are coerced to INTEGER on the way in. A text/integer mismatch on
    either side of the join matches nothing in SQLite while raising no error,
    which would look exactly like a catalogue that simply had no overlap.
    """
    label = label or table
    existing = _row_count(con, table)
    if existing and not reload_:
        print(f"{table}: {existing:,} rows already loaded -- skipping "
              f"(pass --reload to rebuild)")
        return existing

    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"tarva analysispipe: {label} file not found: {path}")

    print(f"Loading {path.name} into {table} (chunks of {chunk_size:,} rows)")
    db = con.cursor()
    db.execute(f"DROP TABLE IF EXISTS {table}")
    con.commit()

    reader = pd.read_csv(
            path, sep="\t", dtype=str, chunksize=chunk_size, usecols=columns,
            low_memory=False,
            )

    total, dropped = 0, 0
    for n, chunk in enumerate(reader, start=1):
        chunk.columns = [str(c).lstrip("#").strip() for c in chunk.columns]
        missing = [c for c in (chrom_col, pos_col) if c not in chunk.columns]
        if missing:
            raise SystemExit(
                    f"tarva analysispipe: {path} has no column(s) {', '.join(missing)}. "
                    f"Columns found: {', '.join(chunk.columns)}"
                    )

        chunk[chrom_col] = chunk[chrom_col].astype(str).str.strip()
        positions = pd.to_numeric(chunk[pos_col], errors="coerce")
        bad = int(positions.isna().sum())
        if bad:
            chunk = chunk.loc[positions.notna()].copy()
            positions = positions.loc[positions.notna()]
            dropped += bad
        chunk[pos_col] = positions.astype("int64")

        chunk.to_sql(table, con, index=False, if_exists="append")
        total += len(chunk)
        print(f"  {table}: chunk {n} -- {total:,} rows loaded", flush=True)

    con.commit()
    db.execute(f"CREATE INDEX IF NOT EXISTS {table}_pos_idx ON {table}({chrom_col},{pos_col});")
    con.commit()
    print(f"{table}: {total:,} rows, indexed on ({chrom_col}, {pos_col})")
    if dropped:
        print(f"WARNING: {dropped:,} row(s) of {path.name} had a non-numeric {pos_col} "
              f"and were skipped", file=sys.stderr)
    return total


#----------------------------------------------------
# the check itself
#----------------------------------------------------

def _attach_tarva(con, db_path):
    db_path = Path(db_path)
    if not db_path.is_file():
        raise SystemExit(f"tarva analysispipe: database not found: {db_path}")
    # read-only: this step must not be able to alter the source of truth, and a
    # read-only attachment also leaves no -wal/-shm beside TARVa.db.
    con.execute("ATTACH DATABASE ? AS tarva", (db_path.resolve().as_uri() + "?mode=ro",))


def _tarva_columns(con, table):
    db = con.cursor()
    present = [c[1] for c in db.execute(f"PRAGMA tarva.table_info({table})").fetchall()]
    if not present:
        raise SystemExit(f"tarva analysispipe: table {table} not found in the database")
    return present


def require_columns(con, table, chrom_col, pos_col):
    present = _tarva_columns(con, table)
    missing = [c for c in (chrom_col, pos_col, "a1type", "rids", "con") if c not in present]
    if missing:
        raise SystemExit(
                f"tarva analysispipe: {table} has no column(s) {', '.join(missing)}. "
                f"Columns found: {', '.join(present)}. Set --chrom-column / --pos-column "
                f"if this table names them differently."
                )


def count_ag_positions(con, table, a1type, chrom_col, pos_col):
    """
    How many distinct (chrom, pos) pairs carry an `a1type` call at all.

    Taken from the database alone, before either catalogue is consulted, so it
    stands as the denominator the match counts are read against rather than
    being a by-product of the join.
    """
    db = con.cursor()
    count_ag = f"""
        SELECT COUNT(DISTINCT t.{chrom_col} || ':' || CAST(t.{pos_col} AS INTEGER))
        FROM tarva.{table} t
        WHERE t.a1type = ?
    """
    return db.execute(count_ag, (a1type,)).fetchone()[0]


def check_positions(con, table, a1type, chrom_col, pos_col, chunk_size,
                    csv_path, result_table=RESULT_TABLE):
    """
    Every distinct (chrom, pos) carrying an `a1type` call, with the number of
    samples per condition holding a call there, and whether each catalogue
    knows the position.

    Both catalogues are tested for every position rather than the second only
    where the first missed, because the overlap between them is itself part of
    the answer.

    The per-condition sample counts are COUNT(DISTINCT rids) per `con` value,
    pivoted into one column per condition inside SQL. That keeps a position to
    exactly one output row, so reading the result in chunks cannot split a
    position's conditions across two chunks.
    """
    require_columns(con, table, chrom_col, pos_col)
    db = con.cursor()

    get_cons = f"""
        SELECT DISTINCT con
        FROM tarva.{table}
        WHERE con IS NOT NULL
        ORDER BY con
    """
    conditions = [c[0] for c in db.execute(get_cons).fetchall()]
    if not conditions:
        raise SystemExit(
                f"tarva analysispipe: no non-null con values in {table} -- nothing to "
                f"count per condition"
                )
    print("Conditions found: {}".format(", ".join(map(str, conditions))))

    # The condition alias stays quoted: its text comes from the data, so an
    # unquoted alias would break on any label that is not a bare identifier.
    count_columns, params = [], []
    for cond in conditions:
        count_columns.append(
                f'COUNT(DISTINCT CASE WHEN t.con = ? THEN t.rids END) AS "n_{cond}"')
        params.append(cond)
    params.append(a1type)
    counts = ",\n                   ".join(count_columns)

    get_checks = f"""
        WITH ag AS (
            SELECT t.{chrom_col} AS chrom, CAST(t.{pos_col} AS INTEGER) AS pos,
                   {counts}
            FROM tarva.{table} t
            WHERE t.a1type = ?
            GROUP BY t.{chrom_col}, CAST(t.{pos_col} AS INTEGER)
        )
        SELECT ag.*,
               EXISTS (SELECT 1 FROM {REDI_TABLE} r
                        WHERE r.{REDI_CHROM} = ag.chrom AND r.{REDI_POS} = ag.pos) AS redi_hit,
               (SELECT GROUP_CONCAT(DISTINCT e.{EDITOME_ID}) FROM {EDITOME_TABLE} e
                 WHERE e.{EDITOME_CHROM} = ag.chrom AND e.{EDITOME_POS} = ag.pos) AS editome_ids
        FROM ag
        ORDER BY ag.chrom, ag.pos
    """

    count_names = [f"n_{c}" for c in conditions]
    out_columns = ["chrom", "pos"] + count_names + ["redi", "editome", "ADeditom_ID"]

    totals = {"positions": 0, "redi": 0, "editome": 0, "both": 0, "neither": 0}
    first = True
    for chunk in pd.read_sql_query(get_checks, con, params=params, chunksize=chunk_size):
        ids = chunk["editome_ids"].fillna("").astype(str).str.strip()
        chunk["redi"] = chunk["redi_hit"].fillna(0).astype(int).map({1: "yes", 0: "no"})
        chunk["editome"] = (ids != "").map({True: "yes", False: "no"})
        # NA, not blank: an empty cell reads as a missing value in R and pandas
        # alike, and here it is a definite negative rather than an absent answer.
        chunk["ADeditom_ID"] = ids.where(ids != "", "NA")
        chunk = chunk[out_columns]

        totals["positions"] += len(chunk)
        totals["redi"] += int((chunk["redi"] == "yes").sum())
        totals["editome"] += int((chunk["editome"] == "yes").sum())
        totals["both"] += int(((chunk["redi"] == "yes") & (chunk["editome"] == "yes")).sum())
        totals["neither"] += int(((chunk["redi"] == "no") & (chunk["editome"] == "no")).sum())

        chunk.to_csv(csv_path, mode="w" if first else "a", header=first, index=False)
        chunk.to_sql(result_table, con, index=False,
                     if_exists="replace" if first else "append")
        first = False
        print(f"  checked {totals['positions']:,} position(s)", flush=True)

    if first:
        empty = pd.DataFrame(columns=out_columns)
        empty.to_csv(csv_path, index=False)
        empty.to_sql(result_table, con, index=False, if_exists="replace")
        print(f"WARNING: no {a1type} calls found in {table} -- nothing to check",
              file=sys.stderr)
    con.commit()
    return totals, conditions


#----------------------------------------------------
# CLI
#----------------------------------------------------

def add_arguments(parser):
    parser.add_argument("--db", required=True, help="TARVa SQLite database (read-only here)")
    parser.add_argument(
            "--redi", required=True, metavar="PATH",
            help=(f"Tab-delimited normal-editing catalogue (REDIportal). Only its distinct "
                  f"{REDI_CHROM}/{REDI_POS} pairs are stored."),
            )
    parser.add_argument(
            "--editome", required=True, metavar="PATH",
            help=(f"Tab-delimited disease-editing catalogue (AD editome), joined on its "
                  f"{EDITOME_CHROM}/{EDITOME_POS} columns"),
            )
    parser.add_argument("--output", required=True, help="Directory to write results into")
    parser.add_argument(
            "--external-db", default=None, metavar="PATH",
            help=f"Where to build the external-source database (default: <output>/{DEFAULT_EXTERNAL_DB})",
            )
    parser.add_argument(
            "--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE,
            help=f"Rows per chunk when reading the source files and the result "
                 f"(default: {DEFAULT_CHUNK_SIZE:,})",
            )
    parser.add_argument(
            "--editome-columns", default=None, metavar="LIST",
            help="Comma-separated columns to keep from --editome (default: all)",
            )
    parser.add_argument(
            "--reload", action="store_true",
            help="Rebuild both source tables even if they are already populated",
            )
    parser.add_argument(
            "--table", default=tarva_db.DEFAULT_TABLE,
            help=f"TARVa table holding the calls (default: {tarva_db.DEFAULT_TABLE})",
            )
    parser.add_argument(
            "--a1type", default=DEFAULT_A1TYPE,
            help=(f"a1type label to check (default: {DEFAULT_A1TYPE} -- the build's label "
                  f"for an A->G substitution, strand already resolved)"),
            )
    parser.add_argument(
            "--chrom-column", default=DEFAULT_CHROM_COLUMN,
            help=f"Chromosome column of --table (default: {DEFAULT_CHROM_COLUMN})",
            )
    parser.add_argument(
            "--pos-column", default=DEFAULT_POS_COLUMN,
            help=f"Position column of --table (default: {DEFAULT_POS_COLUMN})",
            )
    parser.add_argument(
            "--hpc", action="store_true",
            help="Submit this step as a SLURM job instead of running locally",
            )
    parser.add_argument(
            "--conda-env", default="tarva",
            help="Conda environment to activate in the generated SLURM script (default: tarva)",
            )
    parser.add_argument(
            "--directives-dict", dest="directives_dict", type=parse_directives_dict, default={},
            metavar="DICT",
            help=("SLURM sbatch directives as a Python dict literal, only used with --hpc. "
                  "See `tarva qc --help` for recognized keys."),
            )


def run(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.hpc:
        _submit_hpc(args)
        return

    external_db = Path(args.external_db) if args.external_db else output_dir / DEFAULT_EXTERNAL_DB
    external_db.parent.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "09_ag_external_check.csv"

    con = connect_external(external_db)
    try:
        redi_rows = load_redi(con, args.redi, args.chunk_size, reload_=args.reload)
        editome_rows = load_source(
                con, args.editome, EDITOME_TABLE, EDITOME_CHROM, EDITOME_POS, args.chunk_size,
                columns=_csv_list(args.editome_columns), reload_=args.reload, label="--editome")

        _attach_tarva(con, args.db)
        require_columns(con, args.table, args.chrom_column, args.pos_column)
        ag_positions = count_ag_positions(
                con, args.table, args.a1type, args.chrom_column, args.pos_column)
        print(f"distinct {args.a1type} position(s) in {args.table}: {ag_positions:,}")

        totals, conditions = check_positions(
                con, args.table, args.a1type, args.chrom_column, args.pos_column,
                args.chunk_size, csv_path)
    finally:
        con.close()

    summary = {
            "db": str(Path(args.db).resolve()),
            "table": args.table,
            "a1type": args.a1type,
            "external_db": str(external_db.resolve()),
            "redi_file": str(Path(args.redi).resolve()),
            "redi_unique_positions": redi_rows,
            "editome_file": str(Path(args.editome).resolve()),
            "editome_rows": editome_rows,
            "conditions": [str(c) for c in conditions],
            "ag_positions_total": ag_positions,
            "positions_checked": totals["positions"],
            "in_redi": totals["redi"],
            "in_editome": totals["editome"],
            "in_both": totals["both"],
            "in_neither": totals["neither"],
            }
    summary_path = output_dir / "09_external_check_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(f"\ndistinct {args.a1type} position(s) in {args.table}: {ag_positions:,}")
    print(f"  checked:                  {totals['positions']:,}")
    print(f"  in REDIportal (normal):   {totals['redi']:,}")
    print(f"  in the AD editome:        {totals['editome']:,}")
    print(f"  in both:                  {totals['both']:,}")
    print(f"  in neither:               {totals['neither']:,}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {summary_path}")


def _submit_hpc(args):
    output_dir = Path(args.output)
    command_parts = [
            # -u: unbuffered, so a job killed at the wall clock still leaves
            # its progress in the .o rather than discarding the buffer.
            "python", "-u", "-m", "tarva.cli", "analysispipe", "external-checks",
            "--db", str(args.db),
            "--redi", str(args.redi),
            "--editome", str(args.editome),
            "--output", str(args.output),
            "--chunk-size", str(args.chunk_size),
            "--table", str(args.table),
            "--a1type", str(args.a1type),
            "--chrom-column", str(args.chrom_column),
            "--pos-column", str(args.pos_column),
            ]
    if args.external_db:
        command_parts += ["--external-db", str(args.external_db)]
    if args.editome_columns:
        command_parts += ["--editome-columns", str(args.editome_columns)]
    if args.reload:
        command_parts.append("--reload")

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    script_path = write_sbatch_script(
            job_name="tarva_external_checks",
            commands=" ".join(command_parts),
            log_dir=output_dir / "logs",
            script_path=output_dir / "tarva_external_checks.slurm",
            modules=["anaconda3"],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
