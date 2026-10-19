"""
dbpipe load-expression: load per-transcript FPKM and TPM into the expression
table.

TPM is computed over every row of a sample's table before any filtering, since
it is a within-sample composition (TPM_i = FPKM_i / sum(FPKM) * 1e6) and
filtering first would renormalise it over a subset.

    `tarva dbpipe load-expression --stringtie-dir <assembly_dir> --db <TARVa.db>
        [--overwrite]`

Input:  one Stringtie table per individual at
        <stringtie-dir>/<condition>/<rid>/t_data.ctab, as `stringtiepipe
        assemble` round 3 writes with -B, and a database that `annotate` and
        `load-samples` have written.
Output: the expression table, keyed per transcript, and load_expression.* keys
        in run_meta.
"""

import gzip

from pathlib import Path

import pandas as pd

from tarva.dbpipe import schema
from tarva.dbpipe.annotate import split_versions

LOAD_EXPRESSION_TABLES = ("expression",)

CTAB_NAME = "t_data.ctab"
TRANSCRIPT_COLUMN = "t_name"
FPKM_COLUMN = "FPKM"
WRITE_CHUNK = 100000


def add_arguments(parser):
    parser.add_argument("--stringtie-dir", required=True, metavar="DIR",
                        help="The --output directory of `stringtiepipe assemble`. Each "
                             "individual's table is read from "
                             f"<dir>/<condition>/<rid>/{CTAB_NAME}, with a recursive "
                             f"search as a fallback")
    parser.add_argument("--db", required=True,
                        help="Database that `annotate` and `load-samples` have written")
    parser.add_argument("--overwrite", action="store_true",
                        help="Drop and reload expression instead of failing")


#----------------------------------------------------
# inputs
#----------------------------------------------------

def _open_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def load_cohort(con):
    """
    [(individualID, condition)] from the individual table.

    Column names come from run_meta.
    """
    table = schema.INDIVIDUAL_TABLE
    schema.require_table(con, table, "tarva dbpipe load-samples")
    id_column, condition_column = schema.individual_columns(con)
    query = (f'SELECT "{id_column}", "{condition_column}" FROM {table} '
             f'ORDER BY "{id_column}"')
    return [(str(individual), str(condition))
            for individual, condition in con.execute(query)]


def load_known_transcripts(con):
    schema.require_table(con, "transcript", "tarva dbpipe annotate")
    return {str(row[0]) for row in con.execute("SELECT tid FROM transcript")}


def resolve_ctab(root, rid, condition):
    """The documented path first, then a recursive search within that rid."""
    expected = Path(root) / condition / rid / CTAB_NAME
    if expected.is_file():
        return expected
    found = sorted(p for p in Path(root).rglob(CTAB_NAME) if rid in p.parts)
    if len(found) == 1:
        schema.warn(f"{rid}: not at {expected}, using {found[0]}")
        return found[0]
    if len(found) > 1:
        raise SystemExit(
                f"tarva dbpipe load-expression: {len(found)} {CTAB_NAME} files under "
                f"{root} for {rid} ({', '.join(str(p) for p in found[:3])}). Point "
                f"--stringtie-dir at one run's output directory.")
    raise SystemExit(
            f"tarva dbpipe load-expression: no {CTAB_NAME} for {rid} under {root}. "
            f"Expected {expected}. Did `stringtiepipe assemble` run round 3 with -B?")


#----------------------------------------------------
# one sample
#----------------------------------------------------

def read_ctab(path, rid, known_tids):
    """
    (frame, counts) for one individual: rid, tid, FPKM, TPM.

    TPM is computed before the transcript filter, over every quantified row,
    because it is a within-sample composition rather than a per-row value.
    """
    with _open_text(path) as handle:
        table = pd.read_csv(handle, sep="\t", dtype=str)
    table.columns = [str(c).strip() for c in table.columns]

    for column in (TRANSCRIPT_COLUMN, FPKM_COLUMN):
        if column not in table.columns:
            raise SystemExit(
                    f"tarva dbpipe load-expression: {path} has no {column!r} column. "
                    f"Columns found: {', '.join(table.columns)}. Is this a Ballgown "
                    f"{CTAB_NAME}?")

    stems, _ = split_versions(table[TRANSCRIPT_COLUMN].astype(str))
    frame = pd.DataFrame({
        "tid": stems,
        "FPKM": pd.to_numeric(table[FPKM_COLUMN], errors="coerce"),
        })
    frame["individualID"] = rid

    counts = dict(rows=len(frame), no_fpkm=0, unknown_tid=0, written=0)

    no_fpkm = frame["FPKM"].isna()
    counts["no_fpkm"] = int(no_fpkm.sum())
    frame = frame.loc[~no_fpkm].copy()

    # The TPM denominator: the whole quantified transcriptome for this sample.
    total = frame["FPKM"].sum()
    frame["TPM"] = (frame["FPKM"] / total * 1e6) if total > 0 else 0.0
    if total <= 0:
        schema.warn(f"{rid}: total FPKM is {total} -- every TPM is 0. Did round 3 "
                    f"quantify anything?")

    before = len(frame)
    frame = frame.loc[frame["tid"].isin(known_tids)].copy()
    counts["unknown_tid"] = before - len(frame)

    duplicates = frame["tid"].duplicated()
    if int(duplicates.sum()):
        schema.warn(f"{rid}: {int(duplicates.sum()):,} duplicate transcript id(s) in "
                    f"{path.name} -- keeping the first of each")
        frame = frame.loc[~duplicates].copy()

    counts["written"] = len(frame)
    return frame, counts


#----------------------------------------------------
# writing
#----------------------------------------------------

def check_reloadable(con, overwrite):
    existing = con.execute("SELECT COUNT(*) FROM expression").fetchone()[0]
    if existing and not overwrite:
        raise SystemExit(
                f"tarva dbpipe load-expression: expression already holds {existing:,} "
                f"row(s). Pass --overwrite to reload it.")


def prepare_tables(con):
    existing = con.execute("SELECT COUNT(*) FROM expression").fetchone()[0]
    if existing:
        # DROP and recreate, as the other steps do, so --overwrite survives a
        # schema change rather than appending into a stale definition.
        con.execute("DROP TABLE IF EXISTS expression")
        con.commit()
        schema.create_schema(con, LOAD_EXPRESSION_TABLES)
        schema.report("--overwrite: dropped and recreated expression")


#----------------------------------------------------
# entry point
#----------------------------------------------------

def run(args):
    db_path = Path(args.db)
    if not db_path.is_file():
        raise SystemExit(
                f"tarva dbpipe load-expression: database not found: {db_path}. Run "
                f"`dbpipe annotate` and `dbpipe load-samples` first.")

    con = schema.connect_for_build(db_path)
    try:
        schema.create_schema(con)
        check_reloadable(con, args.overwrite)
        cohort = load_cohort(con)
        known_tids = load_known_transcripts(con)
        schema.report(f"cohort: {len(cohort)} individual(s) from `individual`")
        schema.report(f"transcripts: {len(known_tids):,} in `transcript`")

        paths = [(rid, resolve_ctab(args.stringtie_dir, rid, condition))
                 for rid, condition in cohort]
        prepare_tables(con)

        totals = dict(rows=0, no_fpkm=0, unknown_tid=0, written=0)
        for rid, path in paths:
            frame, counts = read_ctab(path, rid, known_tids)
            frame.reindex(columns=["individualID", "tid", "FPKM", "TPM"]).to_sql(
                    "expression", con, if_exists="append", index=False,
                    chunksize=WRITE_CHUNK)
            con.commit()
            for key in totals:
                totals[key] += counts[key]
            schema.report(f"  {rid}: {counts['rows']:,} row(s), "
                          f"{counts['written']:,} written, "
                          f"{counts['unknown_tid']:,} transcript(s) not in the "
                          f"annotation, {counts['no_fpkm']:,} without FPKM")

        if not totals["written"]:
            raise SystemExit(
                    f"tarva dbpipe load-expression: nothing loaded. Every transcript in "
                    f"the ctab tables was absent from `transcript` -- check that the "
                    f"assembly these tables quantify is the one annotate read.")

        schema.create_indexes(con, LOAD_EXPRESSION_TABLES)
        schema.set_meta(con, "load_expression", {
            "load_expression.stringtie_dir": Path(args.stringtie_dir).resolve(),
            "load_expression.n_samples": len(cohort),
            "load_expression.n_rows": totals["written"],
            "load_expression.n_ctab_rows": totals["rows"],
            "load_expression.n_dropped_unknown_tid": totals["unknown_tid"],
            "load_expression.n_dropped_no_fpkm": totals["no_fpkm"],
            })
    finally:
        con.close()

    schema.report("")
    schema.report("load-expression complete")
    schema.report(f"  rows                      {totals['written']:,}")
    schema.report(f"  per sample                "
                  f"{totals['written'] // max(len(cohort), 1):,}")
    schema.report(f"  dropped, not annotated    {totals['unknown_tid']:,}")
    schema.report(f"  dropped, no FPKM          {totals['no_fpkm']:,}")
    schema.report(f"Wrote {db_path}")
    return 0
