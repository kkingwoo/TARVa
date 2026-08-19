import pandas as pd

def split_genes_by_condition_enrichment(
    input_csv: str = "within_genes_set.csv",
    p_cutoff: float = 0.075
):
    # 1. Read input CSV
    df = pd.read_csv(input_csv)

    # 2. Filter for genes with raw p-value < 0.05
    df_sig = df[df["pvalue"] < p_cutoff].copy()

    # 3. Sum counts and calculate condition proportions
    df_sig["total_counts"] = df_sig["ad_counts"] + df_sig["control_counts"]
    df_sig["ad_prop"] = df_sig["ad_counts"] / df_sig["total_counts"]
    df_sig["control_prop"] = df_sig["control_counts"] / df_sig["total_counts"]

    # 4. Assign enrichment groups
    # > 75% in AD -> AD-enriched
    # > 75% in CN -> CN-enriched
    # 29% - 75% in both -> Shared
    def assign_group(row):
        if row["ad_prop"] > 0.75:
            return "AD-enriched"
        elif row["control_prop"] > 0.75:
            return "CN-enriched"
        else:
            return "Shared"

    df_sig["group"] = df_sig.apply(assign_group, axis=1)

    # 5. Format output column headers
    df_out = df_sig.rename(columns={
        "gene": "ENSGID",
        "pvalue": "pval",
        "adj_pvalue": "pval_adj",
        "ad_counts": "ad counts",
        "control_counts": "control counts"
    })

    output_cols = ["ENSGID", "pval", "pval_adj", "ad counts", "control counts"]

    # 6. Export to 3 distinct CSV files
    ad_df = df_out[df_out["group"] == "AD-enriched"][output_cols]
    cn_df = df_out[df_out["group"] == "CN-enriched"][output_cols]
    shared_df = df_out[df_out["group"] == "Shared"][output_cols]

    ad_df.to_csv("ad_enriched_genes.csv", index=False)
    cn_df.to_csv("cn_enriched_genes.csv", index=False)
    shared_df.to_csv("shared_genes.csv", index=False)

    print(f"Export Complete:")
    print(f" - ad_enriched_genes.csv: {len(ad_df)} genes")
    print(f" - cn_enriched_genes.csv: {len(cn_df)} genes")
    print(f" - shared_genes.csv: {len(shared_df)} genes")


if __name__ == "__main__":
    split_genes_by_condition_enrichment()

    df = pd.read_csv("ad_enriched_genes.csv")

    # Filter for rows where the ID starts with 'ENSG'
    ensg_ids = df[df["ENSGID"].astype(str).str.startswith("ENSG")]["ENSGID"].tolist()

    # Print list of ENSG IDs
    for ensg in ensg_ids:
        print(ensg.split('_')[0])
