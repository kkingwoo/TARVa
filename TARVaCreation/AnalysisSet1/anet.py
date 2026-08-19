"""
anet.py

End-to-end RNA-variant stoichiometry pipeline, reading directly from
`sample_tab` (the raw per-position variant call table) in the TARVa
sqlite database.


Pipeline steps, run in order:
  1. within_gene_test   - per gene, ALL call-types pooled, UNWEIGHTED
                           proportions (a1_ad/dp). Sum per sample first,
                           then Mann-Whitney U (AD sums vs Control sums),
                           then Benjamini-Hochberg FDR correction, plus a
                           rank-biserial correlation (rbc) effect size.
  2. classify_and_annotate - keep genes with raw pvalue <= 0.05 (same
                           filter classify_genes.py applied), classify each
                           into AD-enriched / CN-enriched / Shared from its
                           sample_tab site distribution, and look up a HUGO
                           gene name from All_Info_TopGenes_tab (falling
                           back to 'Unknown'). Saved as within_genes_set.csv
                           -- this is the "within genes" csv, now carrying
                           the gene group and gene name columns.
  3. global_calltype_test - per call-type, pooled across the significant
                           genes from step 2, WEIGHTED proportions
                           (a1_ad/dp * 1/gene_len), also with an rbc
                           effect size column.
  4. export_sample_site_counts - per-sample/per-call-type site counts for
                           the same significant gene set.

NOTE on cohort sizes: n_ad/n_control default to COUNT(DISTINCT rid) from
sample_tab, which will UNDERCOUNT if any subject has zero rows in
sample_tab (e.g. zero modified sites anywhere). Since the true cohort is
29 AD / 30 Control, N_AD/N_CONTROL are passed explicitly below until an
authoritative subject-count source (e.g. the clinical metadata CSV) is
wired in.
"""
from __future__ import annotations

import os
import sys

import sqlite3 as sql
import pandas as pd

from scipy.stats import mannwhitneyu
from concurrent import futures

from typing import Optional

TABLE = "sample_tab"

# --- edit-type classification (from edit_typing.py, via make_reviewer_csvs_data.py) ---

EDIT_TYPES = [
    "A-G", "A-C", "A-T", "C-T", "C-G", "C-A", "G-A", "G-T", "G-C",
    "T-A", "T-C", "T-G", "ins", "del",
]
_NEG_STRAND_CALL = ("A", "C", "T", "G")
_NEG_STRAND_RNA = ("T", "G", "A", "C")


def classify_edit_type(strand: str, ref: str, alt1: str) -> Optional[str]:
    """
    Return the edit-type label (one of EDIT_TYPES) for a single variant call,
    or None if the call doesn't match any recognized category.
    """
    if alt1 == "*" or len(ref) > len(alt1):
        return "del"
    if len(ref) < len(alt1):
        return "ins"

    if strand == "+":
        
        if len(ref) == len(alt1):
            label = f"{ref}-{alt1}"
            return label if label in EDIT_TYPES else None
        return None

    if strand == "-":

        if len(ref) == len(alt1):
            ri = _NEG_STRAND_CALL.index(ref) if ref in _NEG_STRAND_CALL else -1
            ai = _NEG_STRAND_CALL.index(alt1) if alt1 in _NEG_STRAND_CALL else -1
            if ri == -1 or ai == -1:
                return None
            ro, ao = _NEG_STRAND_RNA[ri], _NEG_STRAND_RNA[ai]
            label = f"{ro}-{ao}"
            return label if label in EDIT_TYPES else None
        return None

    return None


# --- stats helpers ---

def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Standard BH step-up FDR correction. Returns adjusted p-values in the
    same order as the input list."""
    n = len(pvalues)
    if n == 0:
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


def _cohort_sizes(conn: sqlite3.Connection, n_ad: Optional[int], n_control: Optional[int]) -> tuple[int, int]:
    if n_ad is None:
        n_ad = conn.execute("SELECT COUNT(DISTINCT rid) FROM sample_tab WHERE condition='AD'").fetchone()[0]
    if n_control is None:
        n_control = conn.execute("SELECT COUNT(DISTINCT rid) FROM sample_tab WHERE condition='Control'").fetchone()[0]
    return n_ad, n_control


# --- step 1: within-gene stoichiometry test ---

def within_gene_test(
    db_path: str,
    n_ad: Optional[int] = None,
    n_control: Optional[int] = None,
) -> pd.DataFrame:
    """
where {sample} which will be parsed from metadata file to make each [unique] vcf name)    Per-gene test of modification stoichiometry between AD and Control,
    pooling ALL variant call-types together, using UNWEIGHTED proportions
    (a1_ad / dp).

    Summation logic:
      1. Every row in sample_tab with a non-null ref_id and dp > 0
         contributes a raw proportion = a1_ad / dp.
      2. For each (gene, sample) pair, ALL of that sample's raw proportions
         within that gene -- across every position and every call-type --
         are summed into ONE number for that sample.
      3. Samples with zero modified sites in that gene get 0.0
         (zero-padded to the full cohort size).
      4. Mann-Whitney U (two-sided) compares the AD sample-level sums
         against the Control sample-level sums, for that gene, and a
         rank-biserial correlation (rbc) effect size is computed from the
         U statistic.
      5. Benjamini-Hochberg FDR correction across all genes' p-values.

    ad_counts / control_counts = RAW SITE COUNTS: total number of
    modified-position rows contributing to that gene, per condition --
    NOT the number of samples. Tracked separately from the per-sample
    sums used in the actual test.
    """
    conn = sqlite3.connect(db_path)
    n_ad, n_control = _cohort_sizes(conn, n_ad, n_control)

    rows = conn.execute(
        f"SELECT ref_id, rid, condition, a1_ad, dp FROM {TABLE} "
        f"WHERE ref_id IS NOT NULL AND dp IS NOT NULL AND dp > 0 AND a1_ad IS NOT NULL"
    ).fetchall()
    conn.close()

    per_sample_sums: dict = {}
    site_counts: dict = {}
    for ref_id, rid, condition, a1_ad, dp in rows:
        prop = a1_ad / dp
        per_sample_sums.setdefault(ref_id, {"AD": {}, "Control": {}})
        per_sample_sums[ref_id][condition].setdefault(rid, 0.0)
        per_sample_sums[ref_id][condition][rid] += prop

        site_counts.setdefault(ref_id, {"AD": 0, "Control": 0})
        site_counts[ref_id][condition] += 1

    gene_ids, raw_pvalues, ad_counts_list, control_counts_list, effect_size_list = [], [], [], [], []

    for gene_id, by_cond in per_sample_sums.items():
        ad_vals = list(by_cond["AD"].values())
        con_vals = list(by_cond["Control"].values())
        ad_counts = site_counts[gene_id]["AD"]
        con_counts = site_counts[gene_id]["Control"]

        ad_vals = ad_vals + [0.0] * max(0, n_ad - len(ad_vals))
        con_vals = con_vals + [0.0] * max(0, n_control - len(con_vals))

        if not any(ad_vals) and not any(con_vals):
            continue

        u_stat, p_value = mannwhitneyu(ad_vals, con_vals, alternative="two-sided")

        r_effect = 1 - (2 * u_stat) / (len(ad_vals) * len(con_vals))

        gene_ids.append(gene_id)
        raw_pvalues.append(p_value)
        ad_counts_list.append(ad_counts)
        control_counts_list.append(con_counts)
        effect_size_list.append(r_effect)

    adj_pvalues = benjamini_hochberg(raw_pvalues)

    result = pd.DataFrame({
        "gene": gene_ids,
        "pvalue": raw_pvalues,
        "adj_pvalue": adj_pvalues,
        "ad_counts": ad_counts_list,
        "control_counts": control_counts_list,
        "rbc": effect_size_list,
    }).sort_values("pvalue").reset_index(drop=True)

    return result


# --- step 2: classification + gene-name annotation (from classify_genes.py) ---

def process_single_gene(args):
    """
    Worker function run in parallel via ProcessPoolExecutor.
    Calculates the edit distribution group for a single gene from
    type_by_gene.py logic and maps the resulting group into:
      - AD_only / Primarily_AD -> 'AD-enriched'
      - Control_only / Primarily_Control -> 'CN-enriched'
      - Both -> 'Shared'
    """
    db_path, gene_id = args

    conn = sqlite3.connect(db_path, check_same_thread=False)
    db = conn.cursor()
    s_tab = 'sample_tab'

    from_samps = f"""
        SELECT DISTINCT condition, rid, pos, ref_id 
        FROM {s_tab} 
        WHERE ref_id = ? 
        GROUP BY condition, rid, pos, ref_id
    """
    try:
        samps = db.execute(from_samps, (gene_id,)).fetchall()
    except Exception:
        conn.close()
        return gene_id, 'Shared'

    conn.close()

    if not samps:
        return gene_id, 'Shared'

    dist_dict = {}
    for s in samps:
        con, rids, poss = s[0], s[1], str(s[2])
        if poss not in dist_dict:
            dist_dict[poss] = {}
        if con not in dist_dict[poss]:
            dist_dict[poss][con] = []
        if rids not in dist_dict[poss][con]:
            dist_dict[poss][con].append(rids)

    check_edit_dist = []
    check_outs = {
        'AD': {'unique': 0, 'common': 0},
        'Control': {'unique': 0, 'common': 0}
    }

    for positions in dist_dict.keys():
        for conditions in dist_dict[positions].keys():
            if dist_dict[positions][conditions]:
                sampsize = len(dist_dict[positions][conditions])
                check_edit_dist.append((positions, conditions, sampsize))

    for checking in check_edit_dist:
        co, num = checking[1], checking[2]
        if co in check_outs:
            if num == 1:
                check_outs[co]['unique'] += num
            else:
                check_outs[co]['common'] += num

    ad_unique, ad_common = check_outs['AD']['unique'], check_outs['AD']['common']
    control_unique, control_common = check_outs['Control']['unique'], check_outs['Control']['common']

    unique_sums = ad_unique + control_unique
    common_sums = ad_common + control_common

    ad_common_perc, control_common_perc = 0.0, 0.0
    ad_unique_perc, control_unique_perc = 0.0, 0.0

    if common_sums != 0:
        ad_common_perc = float(ad_common / common_sums)
        control_common_perc = float(control_common / common_sums)
    if unique_sums != 0:
        ad_unique_perc = float(ad_unique / unique_sums)
        control_unique_perc = float(control_unique / unique_sums)

    # 5-group classification logic from type_by_gene.py
    if (ad_unique + ad_common) >= 1 and (control_unique + control_common) == 0:
        raw_group = 'AD_only'
    elif (control_unique + control_common) >= 1 and (ad_unique + ad_common) == 0:
        raw_group = 'Control_only'
    elif (ad_unique_perc + ad_common_perc) >= float(0.75):
        raw_group = 'Primarily_AD'
    elif (control_unique_perc + control_common_perc) >= float(0.75):
        raw_group = 'Primarily_Control'
    else:
        raw_group = 'Both'

    # Map to requested 3 groups: AD-enriched, CN-enriched, Shared
    if raw_group in ['AD_only', 'Primarily_AD']:
        final_group = 'AD-enriched'
    elif raw_group in ['Control_only', 'Primarily_Control']:
        final_group = 'CN-enriched'
    else:
        final_group = 'Shared'

    return gene_id, final_group


def classify_and_annotate(
    df: pd.DataFrame,
    db_path: str,
    gene_col: str = 'gene',
    max_workers: int = 32,
) -> pd.DataFrame:
    """
    Takes the within_gene_test output, keeps genes with raw pvalue <= 0.05
    (same filter classify_genes.py applied), and adds:
      - 'gene group' : AD-enriched / CN-enriched / Shared (via
        process_single_gene, run in parallel over sample_tab)
      - 'gene_name'  : HUGO name looked up from All_Info_TopGenes_tab,
        'Unknown' if not found
    """
    df = df[df['pvalue'] <= 0.05].copy()
    df['ensgid'] = df[gene_col].apply(lambda x: str(x).split('_')[0] if pd.notnull(x) else x)

    # Gene name mapping from database, if available
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    gene_name_map = {}
    try:
        name_query = """
            SELECT DISTINCT ref_base, gene_name 
            FROM All_Info_TopGenes_tab 
            WHERE gene_name IS NOT NULL AND gene_name != ''
        """
        for ref_base, gname in c.execute(name_query).fetchall():
            gene_name_map[str(ref_base).split('_')[0]] = gname
            gene_name_map[str(ref_base)] = gname
    except sqlite3.OperationalError:
        pass
    conn.close()

    # Parallel classification across genes using ProcessPoolExecutor
    gene_list = df[gene_col].tolist()
    tasks = [(db_path, g) for g in gene_list]

    group_results = {}
    print(f"Starting parallel processing for {len(gene_list)} genes using {max_workers} workers...", flush=True)
    with futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(process_single_gene, tasks)
        for original_gene_id, group in results:
            group_results[original_gene_id] = group

    df['gene group'] = df[gene_col].map(group_results)
    df['gene_name'] = df['ensgid'].map(gene_name_map).fillna(df[gene_col].map(gene_name_map)).fillna('Unknown')

    return df


# --- step 3: per-call-type test over the significant gene set ---

def global_calltype_test(
    db_path: str,
    gene_ids: list[str],
    n_ad: Optional[int] = None,
    n_control: Optional[int] = None,
    save_csv: Optional[str] = "TARVa_global_variant_calls.csv",
) -> pd.DataFrame:
    """
    Per-call-type test of modification stoichiometry between AD and Control,
    pooled across the given gene_ids (intended to be the significant genes
    from within_gene_test), using WEIGHTED proportions
    (a1_ad/dp * gene_pos_wt).

    Summation logic:
      1. For every row belonging to a gene in gene_ids, classify the call
         (A-G, C-T, ins, del, etc.), and compute the weighted proportion
         a1prop = (a1_ad/dp) * (1/gene_len), where gene_len is derived from
         gtf_tab (gene length == transcript length).
      2. For each (call-type, sample) pair, ALL of that sample's weighted
         proportions of that call-type -- across every position AND EVERY
         GENE in gene_ids -- are summed into ONE number for that sample.
      3. Samples with zero occurrences of that call-type (within the given
         gene set) get 0.0 (zero-padded to full cohort size).
      4. Mann-Whitney U (two-sided) compares AD sample-level sums against
         Control sample-level sums, for that call-type, and a
         rank-biserial correlation (rbc) effect size is computed from the
         U statistic.
      5. Benjamini-Hochberg FDR correction across all tested call-types.
    """
    conn = sqlite3.connect(db_path)
    n_ad, n_control = _cohort_sizes(conn, n_ad, n_control)

    if not gene_ids:
        conn.close()
        return pd.DataFrame(columns=["variant_call_type", "raw_pvalue", "adj_pvalue", "rbc"])

    placeholders = ",".join("?" * len(gene_ids))

    gene_lengths = dict(conn.execute(
        f"SELECT ref_id, MAX(trans_end) - MIN(trans_start) FROM gtf_tab "
        f"WHERE ref_id IN ({placeholders}) GROUP BY ref_id",
        gene_ids,
    ).fetchall())

    rows = conn.execute(
        f"SELECT ref_id, rid, condition, a1_ad, dp, ref, alt1, strand FROM {TABLE} "
        f"WHERE ref_id IN ({placeholders}) AND dp IS NOT NULL AND dp > 5 AND a1_ad IS NOT NULL",
        gene_ids,
    ).fetchall()
    conn.close()

    per_sample_sums: dict = {}
    for ref_id, rid, condition, a1_ad, dp, ref, alt1, strand in rows:
        gene_len = gene_lengths.get(ref_id)
        if not gene_len or gene_len <= 0:
            continue
        a1type = classify_edit_type(strand, ref, alt1)
        if a1type is None:
            continue

        prop = (a1_ad / dp) * (1 / gene_len)
        per_sample_sums.setdefault(a1type, {"AD": {}, "Control": {}})
        per_sample_sums[a1type][condition].setdefault(rid, 0.0)
        per_sample_sums[a1type][condition][rid] += prop

    call_types, raw_pvalues, effect_sizes = [], [], []
    for a1type, by_cond in per_sample_sums.items():
        ad_vals = list(by_cond["AD"].values())
        con_vals = list(by_cond["Control"].values())

        ad_vals = ad_vals + [0.0] * max(0, n_ad - len(ad_vals))
        con_vals = con_vals + [0.0] * max(0, n_control - len(con_vals))

        if not any(ad_vals) and not any(con_vals):
            continue

        u_stat, p_value = mannwhitneyu(ad_vals, con_vals, alternative="two-sided")

        r_effect = 1 - (2 * u_stat) / (len(ad_vals) * len(con_vals))

        call_types.append(a1type)
        raw_pvalues.append(p_value)
        effect_sizes.append(r_effect)

    adj_pvalues = benjamini_hochberg(raw_pvalues)

    result = pd.DataFrame({
        "variant_call_type": call_types,
        "raw_pvalue": raw_pvalues,
        "adj_pvalue": adj_pvalues,
        "rbc": effect_sizes,
    }).sort_values("raw_pvalue").reset_index(drop=True)

    if save_csv:
        result.to_csv(save_csv, index=False)
        print(f"TARVa_global_variant_calls.csv has been created --> {os.path.abspath(save_csv)}")

    return result


# --- step 4: per-sample/per-call-type site counts for the significant gene set ---

def export_sample_site_counts(
    db_path: str,
    gene_ids: list[str],
    save_csv: str = "sample_site_counts_per_variant.csv",
) -> pd.DataFrame:
    """
    Counts the total modified sites per sample across the provided gene_ids for each 
    recognized variant call-type.
    """
    if not gene_ids:
        df_empty = pd.DataFrame(columns=["variant_call", "total sites", "condition"])
        if save_csv:
            df_empty.to_csv(save_csv, index=False)
        return df_empty

    conn = sqlite3.connect(db_path)
    placeholders = ",".join("?" * len(gene_ids))

    rows = conn.execute(
        f"SELECT rid, condition, ref, alt1, strand FROM {TABLE} "
        f"WHERE ref_id IN ({placeholders}) AND dp IS NOT NULL AND dp > 5 AND a1_ad IS NOT NULL",
        gene_ids,
    ).fetchall()
    conn.close()

    # Aggregate site counts per (rid, condition, variant_call_type)
    counts: dict[tuple[str, str, str], int] = {}
    for rid, condition, ref, alt1, strand in rows:
        a1type = classify_edit_type(strand, ref, alt1)
        if a1type is None:
            continue
        key = (rid, condition, a1type)
        counts[key] = counts.get(key, 0) + 1

    records = []
    for (rid, condition, a1type), total_sites in counts.items():
        records.append({
            "variant_call": a1type,
            "total sites": total_sites,
            "condition": condition
        })

    df_out = pd.DataFrame(records)
    if save_csv:
        df_out.to_csv(save_csv, index=False)

    return df_out


if __name__ == "__main__":
    db_path = sys.argv[1]

    # true, known cohort size -- see the cohort-size note above
    N_AD, N_CONTROL = 29, 30

    within = within_gene_test(db_path, n_ad=N_AD, n_control=N_CONTROL)
    print(f"within_gene_test: {len(within)} genes tested")

    annotated = classify_and_annotate(within, db_path, gene_col='gene', max_workers=32)
    annotated.to_csv("within_genes_set.csv", index=False)
    sig_genes = annotated['gene'].tolist()
    print(f"classify_and_annotate: {len(sig_genes)} genes significant at raw p<=0.05, "
          f"classified and annotated -> within_genes_set.csv")

    calltype = global_calltype_test(db_path, sig_genes, n_ad=N_AD, n_control=N_CONTROL)
    print(f"global_calltype_test: {len(calltype)} call-types tested over {len(sig_genes)} significant genes")

    site_counts = export_sample_site_counts(db_path, sig_genes)
    print(f"export_sample_site_counts: {len(site_counts)} rows exported to 'sample_site_counts_per_variant.csv'")
