"""
Standalone utility: how many of the globally-labelled genes have expression?

    python -m tarva.analysispipe.check_gene_groups \\
        --gene-groups <gene_classes.csv> --output /path/to/results_dir

Run it before `confounding`, which compares composition association across the
gene classes: what matters is not the total recovered but whether the classes
are recovered at comparable rates, since uneven loss would confound a
composition difference with a survivorship difference.

Resolution is attempted in two passes: the label file's own gene symbol
against the expression matrix, then -- for symbols that miss -- its ENSG
against the database's ref_base -> gene_name map, in case the two disagree on
naming.

Input:  a CSV of gene class labels, plus 02_gene_tpm_all.csv (the unfiltered
        matrix) and 02_db_gene_map.csv from --output.
Output: 08_gene_group_resolved.csv in --output, which `confounding` consumes.
"""

import argparse
import sys

from pathlib import Path

import pandas as pd


def build_parser():
    p = argparse.ArgumentParser(
            prog="check_gene_groups",
            description="Resolve globally-labelled genes against the expression matrix",
            )
    p.add_argument("--gene-groups", required=True,
                   help="CSV of globally-labelled genes (Table S4)")
    p.add_argument("--output", required=True,
                   help="analysispipe results directory, holding the step 2 outputs")
    p.add_argument("--gene-name-col", default="gene_name",
                   help="Column holding the HGNC symbol (default: gene_name)")
    p.add_argument("--group-col", default="gene_group",
                   help="Column holding the class label (default: gene_group)")
    p.add_argument("--ensg-col", default="ensgid",
                   help="Column holding the unversioned ENSG (default: ensgid)")
    p.add_argument("--unknown", default="Unknown",
                   help="Value marking an unresolved symbol in the label file (default: Unknown)")
    return p


def _read(path, what):
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"check_gene_groups: {what} not found: {path}")
    return pd.read_csv(path)


def main(argv=None):
    args = build_parser().parse_args(argv)
    out = Path(args.output)

    groups = _read(args.gene_groups, "gene-groups file")
    for col in (args.gene_name_col, args.group_col):
        if col not in groups.columns:
            raise SystemExit(
                    f"check_gene_groups: column {col!r} not in {args.gene_groups}. "
                    f"Found: {', '.join(groups.columns)}"
                    )

    tpm = _read(out / "02_gene_tpm_all.csv", "unfiltered expression matrix")
    if "symbol" not in tpm.columns:
        raise SystemExit("check_gene_groups: 02_gene_tpm_all.csv has no 'symbol' column")
    in_matrix = set(tpm["symbol"].astype(str))

    # ENSG -> symbol from the database, used only where the symbol itself misses.
    ensg_to_symbol = {}
    gene_map_path = out / "02_db_gene_map.csv"
    if gene_map_path.is_file():
        gm = pd.read_csv(gene_map_path)
        key = "ref_base_stripped" if "ref_base_stripped" in gm.columns else "ref_base"
        ensg_to_symbol = dict(zip(gm[key].astype(str), gm["gene_name"].astype(str)))
    else:
        print(f"NOTE: {gene_map_path} not found -- skipping the ENSG fallback",
              file=sys.stderr)

    name = groups[args.gene_name_col].astype(str)
    ensg = (groups[args.ensg_col].astype(str)
            if args.ensg_col in groups.columns else pd.Series([""] * len(groups)))

    resolved, status = [], []
    for nm, eg in zip(name, ensg):
        if nm in in_matrix:
            resolved.append(nm)
            status.append("matched")
            continue
        alt = ensg_to_symbol.get(eg)
        if alt and alt in in_matrix:
            resolved.append(alt)
            status.append("matched via ENSG")
            continue
        resolved.append(None)
        status.append("no symbol" if nm == args.unknown else "absent")

    groups = groups.assign(resolved_symbol=resolved, resolution=status)
    groups.to_csv(out / "08_gene_group_resolved.csv", index=False)

    tab = pd.crosstab(groups["resolution"], groups[args.group_col], margins=True)
    print("\n== Resolution by class ==")
    print(tab.to_string())

    usable = groups["resolved_symbol"].notna()
    print("\n== Usable fraction by class ==")
    for grp, sub in groups.groupby(args.group_col):
        n = len(sub)
        k = int(sub["resolved_symbol"].notna().sum())
        print(f"  {grp:<14} {k:>4} / {n:<4} ({100 * k / n:5.1f}%)")
    print(f"  {'TOTAL':<14} {int(usable.sum()):>4} / {len(groups):<4} "
          f"({100 * usable.sum() / len(groups):5.1f}%)")

    rates = groups.groupby(args.group_col)["resolved_symbol"].apply(lambda s: s.notna().mean())
    spread = float(rates.max() - rates.min())
    print(f"\nSpread in recovery rate across classes: {100 * spread:.1f} percentage points")
    if spread > 0.10:
        print("WARNING: classes are recovered at noticeably different rates, so the "
              "surviving genes are not a fair sample of each class. Step 8's "
              "comparison would partly reflect which genes survived.", file=sys.stderr)

    print(f"\nWrote {out / '08_gene_group_resolved.csv'}")


if __name__ == "__main__":
    sys.exit(main())
