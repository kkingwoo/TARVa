"""
analysispipe biotype exclusion: remove genes for what they are, read from an
    annotation table rather than named in the source.

    The table is a TSV with one row per transcript and at least two columns, by
    default `geneID` and `transcriptType`. A gene is excluded when ANY of its
    transcripts carries a type CONTAINING any of the exclusion substrings
    (default: Mt, rRNA, tRNA, misc_RNA). Substrings rather than equality, so
    "Mt" reaches Mt_rRNA and Mt_tRNA and "rRNA" reaches rRNA_pseudogene without
    those spellings being written down anywhere in this pipeline.

    misc_RNA is in the default list because that is where the Pol III
    structural transcripts sit: RN7SL1 and RN7SL2 are typed misc_RNA, not
    rRNA, and a handful of such transcripts can hold a large share of a
    sample's FPKM budget.

    Matching is case-sensitive, as the substrings are given.

    Every distinct type that matched is printed with its gene count. That
    listing is the check on the substrings: a type appearing there that you did
    not intend to remove means a substring reached further than you meant it to.

    Input:  the annotation TSV and the exclusion substrings.
    Output: the excluded gene list, which `gtf-map` writes into the results
            directory and `expression` reads rather than re-deriving.
"""

import os

import pandas as pd

from pathlib import Path


DEFAULT_EXCLUDE_TYPES = "Mt,rRNA,tRNA,misc_RNA,ribozyme"
DEFAULT_EXCLUDE_PREFIXES = "MT-"
DEFAULT_ID_COLUMN = "geneId"
DEFAULT_TYPE_COLUMN = "transcriptType"
DEFAULT_NAME_COLUMN = "transcriptName"
EXCLUSIONS_CSV = "00_excluded_genes.csv"


def csv_list(raw):
    return [v.strip() for v in str(raw).split(",") if v.strip()]


def add_biotype_arguments(parser):
    """Flags for the annotation-driven exclusion."""
    parser.add_argument(
            "--biotype-tsv", default=None, metavar="TSV",
            help=(
                "Tab-separated annotation table whose genes are filtered by transcript "
                "type. Omit it and no gene is excluded, which leaves structural and "
                "mitochondrial RNA in the TPM denominator."
                ),
            )
    parser.add_argument(
            "--biotype-exclude", default=DEFAULT_EXCLUDE_TYPES, metavar="LIST",
            help=(
                f"Comma-separated substrings; a gene is dropped when any of its "
                f"transcript types contains any of them (default: {DEFAULT_EXCLUDE_TYPES}). "
                f"Matched anywhere in the type string and case-sensitively, so 'Mt' "
                f"covers Mt_rRNA and Mt_tRNA. misc_RNA is included because the 7SL RNAs "
                f"are typed that way rather than as rRNA. Every type actually matched "
                f"is printed, which is the check that a substring did not reach too far."
                ),
            )
    parser.add_argument(
            "--biotype-id-column", default=DEFAULT_ID_COLUMN, metavar="NAME",
            help=f"Gene identifier column in --biotype-tsv (default: {DEFAULT_ID_COLUMN})",
            )
    parser.add_argument(
            "--biotype-type-column", default=DEFAULT_TYPE_COLUMN, metavar="NAME",
            help=f"Transcript type column in --biotype-tsv (default: {DEFAULT_TYPE_COLUMN})",
            )
    parser.add_argument(
            "--biotype-name-column", default=DEFAULT_NAME_COLUMN, metavar="NAME",
            help=(
                f"Transcript name column in --biotype-tsv (default: {DEFAULT_NAME_COLUMN}). "
                f"Read only when --biotype-exclude-prefix is non-empty, and required to "
                f"exist in that case."
                ),
            )
    parser.add_argument(
            "--biotype-exclude-prefix", default=DEFAULT_EXCLUDE_PREFIXES, metavar="LIST",
            help=(
                f"Comma-separated prefixes; a gene is dropped when any of its transcript "
                f"names STARTS WITH one of them (default: {DEFAULT_EXCLUDE_PREFIXES}). "
                f"This is the only route to the 13 mitochondrial protein-coding genes, "
                f"which are typed protein_coding and which no type substring can reach. "
                f"The trailing hyphen is load-bearing: 'MT-' is the HGNC prefix for "
                f"mitochondrially encoded genes, whereas a bare 'MT' would also take "
                f"MTOR, MTHFR, MTUS1 and the MT1/MT2 metallothioneins -- nuclear genes, "
                f"astrocyte-enriched, and not what is being removed here. Pass an empty "
                f"string to disable."
                ),
            )

def load_exclusions(tsv_path, patterns, prefixes=(), id_column=DEFAULT_ID_COLUMN,
                    type_column=DEFAULT_TYPE_COLUMN, name_column=DEFAULT_NAME_COLUMN):
    """
    Return (genes, by_type, by_prefix, n_genes_in_table). `genes` is one row
    per excluded gene; `by_type` and `by_prefix` are the audit of what each of
    the two rules actually caught.
    """
    path = Path(tsv_path)
    if not path.is_file():
        raise SystemExit(f"tarva analysispipe: --biotype-tsv not found: {path}")

    # No comment character: these tables often mark the header line with a
    # leading '#', and treating it as a comment would blank out the one row
    # carrying the column names and promote the first gene to a header.
    table = pd.read_csv(path, sep="\t", dtype=str)
    table.columns = [str(c).lstrip("#").strip() for c in table.columns]

    wanted = [id_column, type_column] + ([name_column] if prefixes else [])
    missing = [c for c in wanted if c not in table.columns]
    if missing:
        raise SystemExit(
                f"tarva analysispipe: {path} has no column(s): {', '.join(missing)}. "
                f"Columns found: {', '.join(table.columns)}. Use --biotype-id-column, "
                f"--biotype-type-column and --biotype-name-column if this table names "
                f"them differently."
                )

    # Only these columns are used; every other field in the table is ignored,
    # whatever it holds. A blank name is tolerated rather than dropped, so a
    # row with no transcript name can still be excluded on its type.
    table = table[wanted].dropna(subset=[id_column, type_column])
    table.columns = ["gene_id", "transcript_type"] + (["transcript_name"] if prefixes else [])
    if prefixes:
        table["transcript_name"] = table["transcript_name"].fillna("")
    for col in table.columns:
        table[col] = table[col].str.strip()
    table = table[(table["gene_id"] != "") & (table["transcript_type"] != "")]
    if table.empty:
        raise SystemExit(f"tarva analysispipe: no usable rows in {path}")

    n_all = table["gene_id"].nunique()
    columns = ["gene_id", "gene_id_base", "transcript_types", "matched", "n_transcripts"]
    empty_type = pd.DataFrame(columns=["transcript_type", "n_genes", "n_transcripts"])
    empty_pre = pd.DataFrame(columns=["prefix", "n_genes", "n_transcripts", "names"])

    # Two independent rules: what a transcript is typed as, and what it is
    # named. The second exists solely for the mitochondrial protein-coding
    # genes, which are typed like any other coding gene and so are invisible
    # to the first.
    rules = {p: table["transcript_type"].str.contains(p, regex=False) for p in patterns}
    for p in prefixes:
        rules[f"name:{p}"] = table["transcript_name"].str.startswith(p)

    flags = pd.DataFrame(rules, index=table.index)
    hit = table[flags.any(axis=1)].copy()
    if hit.empty:
        return pd.DataFrame(columns=columns), empty_type, empty_pre, n_all

    labels = list(rules)
    hit["matched"] = flags.loc[hit.index].apply(
            lambda row: "+".join(p for p in labels if row[p]), axis=1)

    genes = (hit.groupby("gene_id")
                .agg(transcript_types=("transcript_type",
                                       lambda s: ",".join(sorted(set(s)))),
                     matched=("matched",
                              lambda s: "+".join(sorted(set("+".join(s).split("+"))))),
                     n_transcripts=("transcript_type", "size"))
                .reset_index())
    # The version suffix is stripped into a second column rather than in place,
    # so a table of versioned identifiers still joins to an unversioned GTF
    # while the identifier as supplied stays on the record.
    genes["gene_id_base"] = genes["gene_id"].str.replace(r"\.\d+$", "", regex=True)

    if patterns:
        typed = hit[flags.loc[hit.index, patterns].any(axis=1)]
        by_type = (typed.groupby("transcript_type")
                        .agg(n_genes=("gene_id", "nunique"),
                             n_transcripts=("gene_id", "size"))
                        .reset_index()
                        .sort_values("n_genes", ascending=False))
    else:
        by_type = empty_type

    # Names are listed, not merely counted: a prefix rule is blunt enough that
    # seeing exactly which transcripts it took is the only real check on it.
    by_prefix = pd.DataFrame([
        {"prefix": p,
         "n_genes": int(hit.loc[flags.loc[hit.index, f"name:{p}"], "gene_id"].nunique()),
         "n_transcripts": int(flags.loc[hit.index, f"name:{p}"].sum()),
         "names": ",".join(sorted(set(
             hit.loc[flags.loc[hit.index, f"name:{p}"], "transcript_name"])))}
        for p in prefixes
    ]) if prefixes else empty_pre

    return genes[columns], by_type, by_prefix, n_all


def resolve_exclusions(args, output_dir):
    """
    Parse --biotype-tsv, write the gene list the R steps read, and report what
    was matched. Returns the path written, or None when no table was supplied.
    """
    if not getattr(args, "biotype_tsv", None):
        print("No --biotype-tsv: no gene is excluded on annotation grounds.")
        return None

    patterns = csv_list(args.biotype_exclude)
    prefixes = csv_list(getattr(args, "biotype_exclude_prefix", ""))
    if not patterns and not prefixes:
        raise SystemExit(
                "tarva analysispipe: --biotype-exclude and --biotype-exclude-prefix are "
                "both empty. Pass at least one, or omit --biotype-tsv to disable the "
                "filter entirely."
                )

    genes, by_type, by_prefix, n_all = load_exclusions(
            args.biotype_tsv, patterns, prefixes=prefixes,
            id_column=args.biotype_id_column,
            type_column=args.biotype_type_column,
            name_column=args.biotype_name_column,
            )

    out = Path(output_dir) / EXCLUSIONS_CSV
    genes.to_csv(out, index=False)

    print(f"\n== Biotype exclusion ({args.biotype_tsv}) ==")
    print(f"type substrings    : {', '.join(patterns) if patterns else '(none)'} "
          f"(case-sensitive, matched anywhere in the type)")
    print(f"name prefixes      : {', '.join(prefixes) if prefixes else '(none)'}")
    print(f"genes excluded     : {len(genes)} of {n_all} in the table")

    if len(by_type):
        print("\ntypes matched:")
        for row in by_type.itertuples():
            print(f"  {row.transcript_type:<28} {row.n_genes:>7} gene(s)  "
                  f"{row.n_transcripts:>8} transcript(s)")
        print("\nEvery type listed above is being removed. One you did not intend")
        print("means a substring reached further than you meant it to.")

    if len(by_prefix):
        print("\nname prefixes matched:")
        for row in by_prefix.itertuples():
            print(f"  {row.prefix:<28} {row.n_genes:>7} gene(s)  "
                  f"{row.n_transcripts:>8} transcript(s)")
            shown = [n for n in row.names.split(",") if n]
            for i in range(0, min(len(shown), 48), 6):
                print("    " + ", ".join(shown[i:i + 6]))
            if len(shown) > 48:
                print(f"    ... and {len(shown) - 48} more")
        print("\nRead the names above. A prefix rule cannot tell a gene family from a")
        print("locus prefix, so this listing is the whole of its safety check.")

    if not len(by_type) and not len(by_prefix):
        print("\nNothing matched. Check the values in the "
              f"{args.biotype_type_column} and {args.biotype_name_column} columns.")
    print(f"\nWrote {out}")
    return out
    

