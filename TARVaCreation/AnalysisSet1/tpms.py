"""
tpms.py

Parses local database to identify variant calls per gene in the global set, then parses stringtie output for each individual.
TPM values are calculated from FPKM values, and the relevant information is parsed for all individuals across both conditions
for each gene to determine correlations with Spearman rank correlaation between TPM values and number of variant calls within 
each condition. Correlation values between conditions are then subject to Fisher transformation in a comparison to determine 
if there is a difference in correlcations between the conditions per gene.

Output are saved in CSV which is used as input for plotting function in TPM...Rmd.

"""

from __future__ import annotations

import argparse
import os
import re

from concurrent.futures import ProcessPoolExecutor, as_completed

import sqlite3 as sql
import numpy as np
import pandas as pd
from scipy.stats import norm,spearmanr

DB_TABLE = "All_Info_TopGenes_tab"
OUT_COLS = [
        "rid","con","ensgid","ensgid_csv","gene_name","gene_group",
        "FPKM","TPM","varcount","novel_locus",
        ]

GTF_ATTR_RE = re.compile(r'(\w+) "([^"]*)"')

# BH-corrections across p-values
def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    if n ==0:
        return []
    indexed = sorted(range(n), key=lambda i: pvalues[i])
    adjusted = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(indexed), start=1):
        original_rank = n - rank + 1
        p = pvalues[i]
        val = min(prev, p * n / original_rank)
        adjusted[i] = val
        prev = val
    return adjusted

# gene id standards
def canonical_id(raw) -> str | None:
    if not isinstance(raw, str) or not raw:
        return None
    toks = raw.replace(".","_").split("_")
    if toks[0].startswith("ENSG"):
        return toks[0]
    if toks[0] == "MSTRG" and len(toks) >= 2:
        return f"MSTRG_{toks[1]}"
    if len(toks) > 1 and toks[-1].isdigit():
        return "_".join(toks[:-1])
    return "_".join(toks)

# GTF -> transcript_id -> canonical ensgid map
def parse_gtf_transcript_gene_map(gtf_path: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = dict(GTF_ATTR_RE.findall(fields[8]))
            t_id = attrs.get("transcript_id")
            ref_gene_id = attrs.get("ref_gene_id")
            if t_id and ref_gene_id:
                mapping[t_id] = canonical_id(ref_gene_id)
    return mapping


# get global gene list
def load_gene_targets(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df.columns = [c.strip() for c in df.columns]
    if "ensgid" not in df.columns:
        raise ValueError(
                f"Expected an 'ensgid' column in {csv_path}, found: {list(df.columns)}"
                )
    keep = [c for c in ("gene","ensgid","gene_group","gene_name") if c in df.columns]
    targets = (
            df[keep]
            .dropna(subset=["ensgid"])
            .drop_duplicates(subset=[c for c in ("gene","ensgid") if c in keep])
            .reset_index(drop=True)
            )
    for optional in ("gene","gene_group","gene_name"):
        if optional not in targets.columns:
            targets[optional] = np.nan

    source_col = "gene" if targets["gene"].notna().any() else "ensgid"
    targets["canon_id"] = targets[source_col].apply(canonical_id)

    return targets[["gene","ensgid","gene_group","gene_name","canon_id"]]

# distinct individuals from local DB
def get_rid_con_pairs(db_path: str, table: str = DB_TABLE) -> list[tuple[str,str]]:
    conn = sql.connect(db_path, check_same_thread=False)
    try:
        rows = conn.execute(f"SELECT DISTINCT rids, con FROM {table}").fetchall()
    finally:
        conn.close()
    return [(rid, con) for rid, con in rows]

# per individual assessment
class GeneQuantifier:

    @staticmethod
    def find_individual_dir(stringtie_dir: str, rid: str) -> str | None:
        candidate = os.path.join(stringtie_dir, rid)
        print(f"[DEBUG] rid={rid!r} -> candidate={candidate!r} | isdir={os.path.isdir(candidate)}", flush=True)
        return candidate if os.path.isdir(candidate) else None

    @staticmethod
    def process_individual(
            rid: str,
            con: str,
            stringtie_dir: str,
            targets_records: list[dict],
            transcript_gene_map: dict[str,str],
            db_path: str,
            table: str = DB_TABLE,
            ) -> tuple[pd.DataFrame, str | None]:
        targets = pd.DataFrame(targets_records)
        empty = pd.DataFrame(columns=OUT_COLS)

        indiv_dir = GeneQuantifier.find_individual_dir(stringtie_dir,rid)
        print(f"[DEBUG] rid={rid!r} con={con!r} -> indiv_dir={indiv_dir!r}", flush=True)
        if indiv_dir is None:
            return empty, f"no directory named '{rid}' found directly under {stringtie_dir}"

        t_data_path = os.path.join(indiv_dir, "t_data.ctab")
        print(f"[DEBUG] rid={rid!r} t_data_path={t_data_path!r} -> exists={os.path.exists(t_data_path)}", flush=True)
        if not os.path.exists(t_data_path):
            return empty, f"found directory {indiv_dir}, but no t_data.ctab inside it"

        t_data = pd.read_csv(
                t_data_path, sep='\t', usecols=["t_name", "gene_id", "gene_name", "FPKM"]
                )
        total_fpkm = t_data["FPKM"].sum()
        t_data["TPM"] = (t_data["FPKM"] / total_fpkm * 1e6) if total_fpkm > 0 else 0.0

        t_data["ensgid_from_gtf"] = t_data["t_name"].map(transcript_gene_map)
        by_ensgid = (
                t_data.dropna(subset=["ensgid_from_gtf"])
                .groupby("ensgid_from_gtf")[["FPKM","TPM"]].sum()
                )

        t_data["canon_locus"] = t_data["gene_id"].apply(canonical_id)
        by_locus = t_data.groupby("canon_locus")[["FPKM","TPM"]].sum()

        merged = targets.copy()
        is_ensg = merged["canon_id"].astype(str).str.startswith("ENSG")

        merged["FPKM"] = np.nan
        merged["TPM"] = np.nan
        merged.loc[is_ensg, "FPKM"] = merged.loc[is_ensg, "canon_id"].map(by_ensgid["FPKM"])
        merged.loc[is_ensg, "TPM"] = merged.loc[is_ensg, "canon_id"].map(by_ensgid["TPM"])
        merged.loc[~is_ensg, "FPKM"] = merged.loc[~is_ensg, "canon_id"].map(by_locus["FPKM"])
        merged.loc[~is_ensg, "TPM"] = merged.loc[~is_ensg, "canon_id"].map(by_locus["TPM"])
        merged["FPKM"] = merged["FPKM"].fillna(0.0)
        merged["TPM"] = merged["TPM"].fillna(0.0)
        merged["novel_locus"] = ~is_ensg

        conn = sql.connect(db_path, check_same_thread=False)
        try:
            q = f"SELECT ref_base, COUNT(*) FROM {table} WHERE rids = ? AND con = ? GROUP BY ref_base"
            rows = conn.execute(q, [rid, con]).fetchall()
        finally:
            conn.close()

        varcount_map: dict[str, float] = {}
        for ref_base, s in rows:
            key = canonical_id(ref_base)
            if key is None:
                continue
            varcount_map[key] = varcount_map.get(key, 0.0) + (s if s is not None else 0.0)

        merged["rid"] = rid
        merged["con"] = con
        merged["varcount"] = merged["canon_id"].map(varcount_map).fillna(0).astype(int)

        merged = merged.rename(columns={"canon_id": "ensgid", "ensgid": "ensgid_csv"})
        return merged[OUT_COLS].reset_index(drop=True), None

# parallel process individuals
def run_pipeline(
    gene_csv: str,
    db_path: str,
    stringtie_dir: str,
    gtf_path: str,
    out_csv: str = "tpm_variant_master.csv",
    max_workers: int = 16,
    flush_threshold: int = 5000,
) -> pd.DataFrame:
    targets = load_gene_targets(gene_csv)
    targets_records = targets.to_dict("records")
    print(f"loaded {len(targets)} target genes from {gene_csv}", flush=True)

    n_novel = targets["canon_id"].astype(str).str.startswith("MSTRG").sum()
    if n_novel:
        print(
            f"[INFO] {n_novel} target genes are StringTie-only loci (no Ensembl "
            f"gene) -- matched by canonical locus id, flagged novel_locus=True, "
            f"and will typically show varcount=0 (see docstring assumption 5).",
            flush=True,
        )

    print(f"parsing GTF: {gtf_path}", flush=True)
    transcript_gene_map = parse_gtf_transcript_gene_map(gtf_path)
    print(f"{len(transcript_gene_map)} transcripts mapped to a reference gene", flush=True)

    pairs = get_rid_con_pairs(db_path)
    print(f"{len(pairs)} individuals found in {DB_TABLE}", flush=True)

    accrued: list[pd.DataFrame] = []
    accrued_rows = 0
    master_chunks: list[pd.DataFrame] = []
    first_write = True

    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = {
            ex.submit(
                GeneQuantifier.process_individual,
                rid, con, stringtie_dir, targets_records,
                transcript_gene_map, db_path,
            ): (rid, con)
            for rid, con in pairs
        }
        for fut in as_completed(futures):
            rid, con = futures[fut]
            try:
                df, reason = fut.result()
            except Exception as e:  
                print(f"[WARN] {rid} ({con}) failed: {e}", flush=True)
                continue
            if df.empty:
                print(f"[WARN] {rid} ({con}): {reason} no matching StringTie dir/data found", flush=True)
                continue

            accrued.append(df)
            accrued_rows += len(df)

            if accrued_rows >= flush_threshold:
                chunk = pd.concat(accrued, ignore_index=True)
                chunk.to_csv(
                    out_csv, mode="w" if first_write else "a",
                    header=first_write, index=False,
                )
                first_write = False
                master_chunks.append(chunk)
                print(f"flushed {len(chunk)} rows -> {out_csv}", flush=True)
                accrued, accrued_rows = [], 0

    if accrued:
        chunk = pd.concat(accrued, ignore_index=True)
        chunk.to_csv(
            out_csv, mode="w" if first_write else "a",
            header=first_write, index=False,
        )
        master_chunks.append(chunk)

    master = (
        pd.concat(master_chunks, ignore_index=True) if master_chunks else pd.DataFrame(columns=OUT_COLS)
    )
    print(f"pipeline complete: {len(master)} total rows -> {out_csv}", flush=True)
    return master

# statistical tests
def tpm_variant_relationship(master: pd.DataFrame) -> pd.DataFrame:
    """
     To answer question: 
        For each gene, across the individuals within a condition, is TPM
        associated with variant count? 
    Then:
      Does that association differ between AD and Control?

    Genes with fewer than 4 individuals, or with no variance in TPM or in
    varcount, in a given condition get rho/p = NaN for that condition
    (Spearman is undefined/unstable there) and are excluded from the FDR
    step rather than silently treated as non-significant.
    """
    records = []
    for ensgid, g in master.groupby("ensgid"):
        gene_name = g["gene_name"].iloc[0] if g["gene_name"].notna().any() else np.nan
        gene_group = g["gene_group"].iloc[0] if g["gene_group"].notna().any() else np.nan
        novel_locus = bool(g["novel_locus"].iloc[0])
        row = {
            "ensgid": ensgid, "gene_name": gene_name, "gene_group": gene_group,
            "novel_locus": novel_locus,
        }

        rs, ns = {}, {}
        for cond in ("AD", "Control"):
            sub = g[g["con"] == cond]
            if len(sub) >= 4 and sub["TPM"].nunique() > 1 and sub["varcount"].nunique() > 1:
                r, p = spearmanr(sub["TPM"], sub["varcount"])
            else:
                r, p = np.nan, np.nan
            rs[cond], ns[cond] = r, len(sub)
            row[f"{cond.lower()}_rho"] = r
            row[f"{cond.lower()}_pvalue"] = p
            row[f"{cond.lower()}_n"] = ns[cond]

        r1, n1 = rs.get("AD", np.nan), ns.get("AD", 0)
        r2, n2 = rs.get("Control", np.nan), ns.get("Control", 0)
        if (
            pd.notna(r1) and pd.notna(r2)
            and n1 > 3 and n2 > 3
            and abs(r1) < 1 and abs(r2) < 1
        ):
            z1, z2 = np.arctanh(r1), np.arctanh(r2)
            se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
            z = (z1 - z2) / se
            p_diff = 2 * (1 - norm.cdf(abs(z)))
        else:
            p_diff = np.nan
        row["diff_pvalue"] = p_diff
        records.append(row)

    result = pd.DataFrame(records)
    result["diff_adj_pvalue"] = np.nan
    valid = result["diff_pvalue"].notna()
    if valid.any():
        result.loc[valid, "diff_adj_pvalue"] = benjamini_hochberg(
            result.loc[valid, "diff_pvalue"].tolist()
        )
    return result.sort_values("diff_pvalue", na_position="last").reset_index(drop=True)


# CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene-csv", required=True, help="target gene list CSV (needs 'ensgid' column)")
    ap.add_argument("--db", required=True, help="path to the sqlite variant-call database")
    ap.add_argument("--stringtie-dir", required=True, help="parent dir of per-individual StringTie output")
    ap.add_argument("--gtf", required=True, help="StringTie guided/merged reference GTF (has ref_gene_id)")
    ap.add_argument("--out-prefix", default="tpm_variant", help="prefix for output CSVs")
    ap.add_argument("--max-workers", type=int, default=16)
    ap.add_argument("--flush-threshold", type=int, default=5000)
    args = ap.parse_args()

    master_csv = f"{args.out_prefix}_master.csv"
    master = run_pipeline(
        gene_csv=args.gene_csv,
        db_path=args.db,
        stringtie_dir=args.stringtie_dir,
        gtf_path=args.gtf,
        out_csv=master_csv,
        max_workers=args.max_workers,
        flush_threshold=args.flush_threshold,
    )

    if master.empty:
        print("no data collected -- check paths/assumptions above before running stats", flush=True)
        return

    stats = tpm_variant_relationship(master)
    stats_csv = f"{args.out_prefix}_relationship_stats.csv"
    stats.to_csv(stats_csv, index=False)
    print(f"relationship stats -> {stats_csv} ({len(stats)} genes tested)", flush=True)


if __name__ == "__main__":
    main()





