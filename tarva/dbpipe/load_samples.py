"""
dbpipe load-samples: load the metadata CSV into the individual table

    `tarva dbpipe load-samples --metadata <meta.csv> --db <TARVa.db>
        [--id-column individualID] [--condition-column group] [--overwrite]`

The whole metadata file becomes the `individual` table, column names and all.
Nothing is renamed, reshaped or dropped.

--id-column and --condition-column only IDENTIFY two columns; they do not
rename them. Both names go into run_meta and are read back through
`schema.individual_columns`. The condition column must hold exactly two
distinct values.

Writes: individual. Records load_samples.* keys in run_meta.
"""

from pathlib import Path

import pandas as pd

from tarva.dbpipe import schema

LOAD_SAMPLES_TABLES = (schema.INDIVIDUAL_TABLE,)
ID_INDEX_NAME = "individual_id_idx"

DEFAULT_ID_COLUMN = "individualID"
DEFAULT_CONDITION_COLUMN = "group"
REQUIRED_CONDITIONS = 2

# What an emptied CSV cell looks like once it has been through astype(str).
BLANK_VALUES = ("", "nan", "NaN", "None", "NA", "na", "<NA>")


def add_arguments(parser):
    parser.add_argument("--metadata", required=True,
                        help="The metadata CSV. Every column becomes a column of the "
                             "individual table, under its own name")
    parser.add_argument("--db", required=True,
                        help="SQLite database to create or add to")
    parser.add_argument("--id-column", default=DEFAULT_ID_COLUMN,
                        help=f"Which column holds the individual identifier (default: "
                             f"{DEFAULT_ID_COLUMN}). The column is NOT renamed. Its "
                             f"values must match the sample names in the VCFs. ")
    parser.add_argument("--condition-column", default=DEFAULT_CONDITION_COLUMN,
                        help=f"Which column holds the comparison variable (default: "
                             f"{DEFAULT_CONDITION_COLUMN}). Not renamed either. Exactly "
                             f"two distinct values are required. ")
    parser.add_argument("--overwrite", action="store_true",
                        help="Replace an existing individual table instead of failing. ")


#----------------------------------------------------
# validation
#----------------------------------------------------

def validate(frame, args):
    """
    Check the three things the rest of the schema relies on, and return the
    two condition labels. Nothing about the frame is modified.
    """
    available = list(frame.columns)
    missing = {flag: name for flag, name in
               (("--id-column", args.id_column),
                ("--condition-column", args.condition_column))
               if name not in available}
    if missing:
        raise SystemExit(
                f"tarva dbpipe load-samples: {args.metadata} has no column for "
                f"{'; '.join(f'{flag} {name!r}' for flag, name in missing.items())}. "
                f"Columns found: {', '.join(available)}")

    if args.id_column == args.condition_column:
        raise SystemExit(
                f"tarva dbpipe load-samples: --id-column and --condition-column are "
                f"both {args.id_column!r}. They name different things.")

    duplicated = [name for name in set(available) if available.count(name) > 1]
    if duplicated:
        raise SystemExit(
                f"tarva dbpipe load-samples: {args.metadata} has duplicate column "
                f"name(s) {', '.join(sorted(duplicated))}. Every column becomes a "
                f"table column, so each must be named once.")

    for flag, column in (("--id-column", args.id_column),
                         ("--condition-column", args.condition_column)):
        values = frame[column].astype(str).str.strip()
        blank = values.isin(BLANK_VALUES)
        if int(blank.sum()):
            raise SystemExit(
                    f"tarva dbpipe load-samples: {int(blank.sum())} row(s) have no "
                    f"value in {flag} ({column!r}). Every individual needs both an "
                    f"identifier and a condition.")

    ids = frame[args.id_column].astype(str).str.strip()
    duplicates = ids.duplicated()
    if int(duplicates.sum()):
        raise SystemExit(
                f"tarva dbpipe load-samples: {int(duplicates.sum())} duplicate "
                f"{args.id_column} value(s) "
                f"(e.g. {', '.join(ids[duplicates].unique()[:3])}). One row per "
                f"individual -- the UNIQUE index created after the write would fail "
                f"anyway.")

    labels = sorted(frame[args.condition_column].astype(str).str.strip().unique())
    if len(labels) != REQUIRED_CONDITIONS:
        raise SystemExit(
                f"tarva dbpipe load-samples: {args.condition_column!r} holds "
                f"{len(labels)} distinct value(s) ({', '.join(labels)}), and exactly "
                f"{REQUIRED_CONDITIONS} are required. Every published test is a "
                f"two-sample Mann-Whitney, so supporting k groups means choosing a "
                f"different test -- new analysis rather than a load option.")
    return labels


#----------------------------------------------------
# writing
#----------------------------------------------------

def write_individual(con, frame, args):
    """
    to_sql with if_exists="replace", then the UNIQUE index on the identifier.

    """
    table = schema.INDIVIDUAL_TABLE
    present = con.execute(
            "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=?",
            (table,)).fetchone()[0]
    if present:
        rows = con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        if rows and not args.overwrite:
            raise SystemExit(
                    f"tarva dbpipe load-samples: {table} already holds {rows:,} row(s). "
                    f"Pass --overwrite to replace it. ")
        if rows:
            schema.report(f"--overwrite: replacing {table} ({rows:,} row(s))")

    frame.to_sql(table, con, if_exists="replace", index=False)
    con.execute(f'DROP INDEX IF EXISTS {ID_INDEX_NAME}')
    con.execute(f'CREATE UNIQUE INDEX {ID_INDEX_NAME} '
                f'ON {table}("{args.id_column}")')
    con.commit()
    schema.report(f"{table}: {len(frame):,} row(s), {len(frame.columns)} column(s) "
                  f"written")
    schema.report(f"  columns: {', '.join(frame.columns)}")
    schema.report(f"  unique index on {args.id_column!r}")
    schema.create_indexes(con, LOAD_SAMPLES_TABLES)


#----------------------------------------------------
# entry point
#----------------------------------------------------

def run(args):
    db_path = Path(args.db)
    metadata_path = Path(args.metadata)
    if not metadata_path.is_file():
        raise SystemExit(
                f"tarva dbpipe load-samples: --metadata not found: {metadata_path}")

    frame = pd.read_csv(metadata_path, dtype=str)
    frame.columns = [str(c).strip() for c in frame.columns]
    schema.report(f"metadata: {len(frame):,} row(s), {len(frame.columns)} column(s) "
                  f"from {metadata_path.name}")

    labels = validate(frame, args)
    sizes = frame.groupby(args.condition_column, sort=True).size()

    con = schema.connect_for_build(db_path)
    try:
        schema.create_schema(con)
        write_individual(con, frame, args)
        schema.set_meta(con, "load_samples", {
            "load_samples.metadata": metadata_path.resolve(),
            "load_samples.id_column": args.id_column,
            "load_samples.condition_column": args.condition_column,
            "load_samples.n_individuals": len(frame),
            "load_samples.columns": ",".join(frame.columns),
            "load_samples.conditions": ",".join(labels),
            "load_samples.condition_sizes": ",".join(
                    f"{label}={count}" for label, count in sizes.items()),
            })
    finally:
        con.close()

    schema.report("")
    schema.report("load-samples complete")
    for label, count in sizes.items():
        schema.report(f"  {label:<26}{count:,}")
    schema.report(f"  total                     {len(frame):,}")
    schema.report(f"Wrote {db_path}")
    return 0
