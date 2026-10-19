"""
Database access for analysispipe.

All SQLite access lives here rather than in the R steps, matching how every
other TARVa table is written, so one module understands the schema.

The sample list comes from All_Info_TopGenes_tab rather than a metadata CSV.
That is deliberate: the cohort is then defined by the same table results are
keyed back to, so the rids join is correct by construction instead of by later
validation. It also needs no FASTQ filename prefixes, which a database built
from an existing call set may not have.

Input:  the database path, and the table names each loader is given.
Output: DataFrames -- the cohort, the gene map, the covariates, the
        per-individual burden and the site-level calls.
"""

import sqlite3
import sys

from pathlib import Path

import pandas as pd

DEFAULT_TABLE = "All_Info_TopGenes_tab"
DEFAULT_CTAB = "t_data.ctab"


def _connect(db_path):
    db_path = Path(db_path)
    if not db_path.is_file():
        raise SystemExit(f"tarva analysispipe: database not found: {db_path}")
    return sqlite3.connect(str(db_path))


def _query(db_path, sql, table, params=()):
    con = _connect(db_path)
    try:
        return con.execute(sql, params).fetchall()
    except sqlite3.OperationalError as exc:
        raise SystemExit(
                f"tarva analysispipe: could not read {table} from {db_path}: {exc}"
                )
    finally:
        con.close()


def load_burden(db_path, table=DEFAULT_TABLE):
    """
    Per-individual editing-call burden.

    The most direct test of the confounding mechanism: if an individual's total
    call count tracks their estimated cell composition, then composition drives
    detection and every count-based finding inherits that. If it does not, the
    mechanism is rebutted at its source rather than at its consequences.
    """
    rows = _query(
            db_path,
            f"SELECT rids, con, COUNT(*) AS n_calls, COUNT(DISTINCT poss) AS n_positions, "
            f"COUNT(DISTINCT ref_base) AS n_genes FROM {table} GROUP BY rids, con",
            table,
            )
    return pd.DataFrame(rows, columns=["rids", "con", "n_calls", "n_positions", "n_genes"])


def load_gene_sites(db_path, genes, table=DEFAULT_TABLE):
    """
    Distinct (gene, individual, position) triples for the given genes -- the
    unit the global enrichment classification counts over, and what the label
    permutation needs.
    """
    genes = [str(g) for g in genes]
    frames = []
    for i in range(0, len(genes), 500):
        chunk = genes[i:i + 500]
        ph = ",".join("?" * len(chunk))
        rows = _query(
                db_path,
                f"SELECT DISTINCT ref_base, rids, poss FROM {table} "
                f"WHERE ref_base IN ({ph})",
                table, tuple(chunk),
                )
        frames.append(pd.DataFrame(rows, columns=["ref_base", "rids", "poss"]))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
            columns=["ref_base", "rids", "poss"])


def load_site_level(db_path, genes, table=DEFAULT_TABLE):
    """
    Per-individual editing proportions at each position of the given genes,
    for the adjusted models. `a1ad`/`dp` is the outcome; `dp` also serves as
    the depth covariate testing the detection pathway directly.
    """
    genes = [str(g) for g in genes]
    frames = []
    for i in range(0, len(genes), 500):
        chunk = genes[i:i + 500]
        ph = ",".join("?" * len(chunk))
        rows = _query(
                db_path,
                f"SELECT ref_base, poss, a1type, rids, con, a1ad, dp FROM {table} "
                f"WHERE ref_base IN ({ph})",
                table, tuple(chunk),
                )
        frames.append(pd.DataFrame(
            rows, columns=["ref_base", "poss", "a1type", "rids", "con", "a1ad", "dp"]))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
            columns=["ref_base", "poss", "a1type", "rids", "con", "a1ad", "dp"])


def load_samples(db_path, table=DEFAULT_TABLE):
    """Sample IDs and group labels -- the cohort definition."""
    rows = _query(db_path, f"SELECT DISTINCT rids, con FROM {table} ORDER BY rids", table)
    if not rows:
        raise SystemExit(f"tarva analysispipe: no samples found in {table} of {db_path}")
    return pd.DataFrame(rows, columns=["rids", "con"])


def load_gene_map(db_path, table=DEFAULT_TABLE):
    """
    ref_base -> gene_name (HUGO). Both the versioned and unversioned forms are
    emitted so either can join; the confounding step needs this to reach
    variant-call data from a symbol.
    """
    rows = _query(
            db_path,
            f"SELECT DISTINCT ref_base, gene_name FROM {table} "
            f"WHERE gene_name IS NOT NULL AND gene_name != ''",
            table,
            )
    df = pd.DataFrame(rows, columns=["ref_base", "gene_name"])
    df["ref_base_stripped"] = (
            df["ref_base"].astype(str).str.replace(r"[._][0-9]+$", "", regex=True)
            )
    return df


DEFAULT_COVARIATES = ("msex", "age_death", "pmi", "braaksc", "ceradsc", "apoe_genotype")
DEFAULT_COVARIATE_KEY = "individualID"


def add_covariate_file_arguments(parser):
    """
    Flags for merging covariates the editing table does not carry. Shared by
    every step that loads covariates, so the two written covariate files stay
    identical to each other.
    """
    parser.add_argument(
            "--covariate-file", default=None, metavar="PATH",
            help=(
                "Delimited file of additional per-sample covariates, merged onto those "
                "read from the database. Needed for variables the editing table does "
                "not hold -- PMI in this cohort. The delimiter is detected, so .csv and "
                ".tsv both work."
                ),
            )
    parser.add_argument(
            "--covariate-file-key", default=DEFAULT_COVARIATE_KEY, metavar="NAME",
            help=(
                f"Column of --covariate-file holding the sample identifier that matches "
                f"rids (default: {DEFAULT_COVARIATE_KEY})"
                ),
            )
    parser.add_argument(
            "--covariate-file-columns", default=None, metavar="LIST",
            help=(
                "Comma-separated columns to take from --covariate-file. Default: every "
                "column except the key. A column the database already supplied is "
                "skipped and reported rather than overwritten."
                ),
            )


def merge_covariate_file(frame, path, key=DEFAULT_COVARIATE_KEY, columns=None):
    """
    Merge extra per-sample covariates from a delimited file onto the frame
    `load_covariates` returned, joining `key` in the file to `rids`.

    This cohort's clinical variables are split across two places: most sit in
    the editing table, but PMI does not, and the file holding it keys on
    individualID. Rather than hardcode that schema, the path, the key column
    and the columns taken are all supplied at runtime.

    The database is treated as authoritative: a column it already provided is
    reported and left alone rather than overwritten from the file.
    """
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"tarva analysispipe: --covariate-file not found: {path}")

    # sep=None with the python engine sniffs the delimiter, so a comma-, tab-
    # or whitespace-separated file all work without being declared. Everything
    # is read as text; the R steps coerce what they need.
    extra = pd.read_csv(path, sep=None, engine="python", dtype=str)
    extra.columns = [str(c).lstrip("#").strip() for c in extra.columns]

    if key not in extra.columns:
        raise SystemExit(
                f"tarva analysispipe: {path} has no column {key!r}. Columns found: "
                f"{', '.join(extra.columns)}. Set --covariate-file-key if the sample "
                f"identifier is named differently."
                )

    wanted = [c for c in (columns if columns else list(extra.columns)) if c != key]
    unknown = [c for c in wanted if c not in extra.columns]
    if unknown:
        raise SystemExit(
                f"tarva analysispipe: {path} has no column(s): {', '.join(unknown)}. "
                f"Columns found: {', '.join(extra.columns)}"
                )

    clash = [c for c in wanted if c in frame.columns]
    if clash:
        print(f"NOTE: {', '.join(clash)} already came from the database -- keeping the "
              f"database value", file=sys.stderr)
        wanted = [c for c in wanted if c not in clash]
    if not wanted:
        return frame

    extra = extra[[key] + wanted].rename(columns={key: "rids"})
    extra["rids"] = extra["rids"].astype(str).str.strip()
    dup = extra["rids"][extra["rids"].duplicated()].unique()
    if len(dup):
        print(
                f"WARNING: {path.name} has repeated rows for {len(dup)} sample(s) "
                f"(e.g. {', '.join(map(str, dup[:3]))}) -- keeping the first",
                file=sys.stderr,
                )
        extra = extra.drop_duplicates("rids")

    frame = frame.copy()
    frame["rids"] = frame["rids"].astype(str)
    merged = frame.merge(extra, on="rids", how="left")

    # Coverage is reported per column rather than assumed. A covariate present
    # for only part of the cohort silently shrinks the adjusted models' sample
    # size, and an identifier-format mismatch matches nothing at all.
    for col in wanted:
        n = int(merged[col].notna().sum())
        print(f"Merged {col} from {path.name}: {n} of {len(merged)} sample(s) matched")
        if n == 0:
            print(
                    f"WARNING: no {key} value in {path.name} matched any rids -- check "
                    f"that the two use the same identifier format",
                    file=sys.stderr,
                    )
    return merged


def _table_columns(db_path, table):
    con = _connect(db_path)
    try:
        return {row[1] for row in con.execute(f"PRAGMA table_info({table})").fetchall()}
    finally:
        con.close()


def load_covariates(db_path, table=DEFAULT_TABLE, columns=DEFAULT_COVARIATES):
    """
    Per-sample clinical covariates. These sit in the same table as the editing
    calls, so no external metadata file is needed -- which matters here, since
    the cohort's FASTQ-era metadata no longer exists.

    Columns are requested individually and a column absent from the table is
    warned about and skipped rather than killing the step, so the default list
    can carry variables that only some databases hold.
    """
    frame = load_samples(db_path, table)
    present = _table_columns(db_path, table)
    for col in columns:
        if col not in present:
            print(f"NOTE: column {col!r} is not in {table} -- skipping", file=sys.stderr)
            continue
        rows = _query(
                db_path,
                f"SELECT DISTINCT rids, {col} FROM {table} WHERE {col} IS NOT NULL",
                table,
                )
        got = pd.DataFrame(rows, columns=["rids", col])
        dup = got["rids"][got["rids"].duplicated()].unique()
        if len(dup):
            print(
                    f"WARNING: {col} is not constant within {len(dup)} sample(s) "
                    f"(e.g. {', '.join(map(str, dup[:3]))}) -- keeping the first value",
                    file=sys.stderr,
                    )
            got = got.drop_duplicates("rids")
        frame = frame.merge(got, on="rids", how="left")
    return frame


def build_crosswalk(stringtie_root, samples, ctab_name=DEFAULT_CTAB):
    """
    Locate <stringtie_root>/<rids>/<ctab_name> for every sample in the cohort.

    Driven by the database, not the filesystem: the cohort defines what must be
    present, so a sample with no directory is fatal (it would silently vanish
    from the matrix and skew composition estimates). Extra directories under
    the root are simply not looked at -- they are quantified samples outside
    the cohort and have no bearing on the analysis.

    Returns (crosswalk_df, missing_rids). Every row of the frame is matched, so
    downstream code never has to filter it.
    """
    root = Path(stringtie_root)
    rows, missing = [], []
    for rid, con in zip(samples["rids"].astype(str), samples["con"]):
        ctab = root / rid / ctab_name
        if ctab.is_file():
            rows.append({"ctab_path": str(ctab), "rids": rid, "con": con})
        else:
            missing.append(rid)
    return pd.DataFrame(rows, columns=["ctab_path", "rids", "con"]), missing
