"""
Figure 5: Forest plot of odds ratios for unique variant calls, by variant call-type
and global gene subset.

    python Fig5_uniqueness_forest.py --input-dir <dir> [--output Fig5.png]

Computes nothing: the odds ratio, its confidence interval and both p-values
come from the step, so the figure and the database cannot disagree.

No main title.
Rows are bolded on whichever p-value column the step actually used, named in
the CSV's `basis` column; when neither basis qualified, nothing is bolded,
because nothing was significant.

Input:  site_uniqueness.csv from `tarva analysispipe site-uniqueness`, in
        --input-dir. Only the `*_enriched` subsets are plotted; the `Both`
        subset is in the CSV and in `uniqueness_test` but is not one of the
        figure's sections.
Output: the image named by --output (the extension picks the format).
"""

import argparse

from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")

import matplotlib.pyplot as plt   # noqa: E402  (after the backend is set)

CSV_NAME = "site_uniqueness.csv"

ALPHA = 0.05

# Okabe-Ito-adjacent: case warm, reference cool.
CASE_COLOUR = "#cc3311"
REF_COLOUR = "#0077bb"

SUBSET_SUFFIX = "_enriched"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True,
                        help="Directory holding site_uniqueness.csv")
    parser.add_argument("--output", default="Fig5_uniqueness_forest.png",
                        help="Image to write (extension picks the format)")
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def load(input_dir):
    """
    The step's output, with the two condition labels and the significance
    basis it chose.
    """
    frame = pd.read_csv(Path(input_dir) / CSV_NAME)
    if frame.empty:
        raise SystemExit(f"{CSV_NAME} is empty -- run `tarva analysispipe "
                         f"site-uniqueness --output {input_dir}` first.")

    for column in ("case_condition", "ref_condition"):
        if column not in frame.columns:
            raise SystemExit(
                    f"{CSV_NAME} has no {column} column. Re-run "
                    f"site-uniqueness: the figure needs the condition labels "
                    f"rather than assuming them.")

    case = str(frame["case_condition"].iloc[0])
    reference = str(frame["ref_condition"].iloc[0])

    basis = [value for value in frame["basis"].dropna().unique() if value]
    if len(basis) > 1:
        raise SystemExit(f"{CSV_NAME} names more than one significance basis: "
                         f"{', '.join(basis)}")
    return frame, case, reference, (basis[0] if basis else None)


def sections(frame, case, reference):
    """
    [(subset, colour, rows)]

    """
    out = []
    for label, colour in ((f"{case}{SUBSET_SUFFIX}", CASE_COLOUR),
                          (f"{reference}{SUBSET_SUFFIX}", REF_COLOUR)):
        rows = frame.loc[frame["edit_dist"] == label]
        rows = rows.loc[rows["odds_ratio"].notna()]
        if rows.empty:
            continue
        out.append((label, colour,
                    rows.sort_values("odds_ratio", ascending=False)))
    return out


def build_sequence(blocks):
    """
    The plot's rows, top to bottom: a header, then its call types, then a
    spacer before the next section.

    """
    sequence = []
    for index, (label, colour, rows) in enumerate(blocks):
        if index:
            sequence.append(("spacer", None, None))
        sequence.append(("header", label, colour))
        for row in rows.itertuples():
            sequence.append(("row", row, colour))
    # Reversed so the first section appears at the top of the axes.
    return sequence[::-1]


def render(sequence, basis, case, reference, output, dpi):
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

    height = max(4.0, 0.32 * len(sequence) + 1.5)
    figure = plt.figure(figsize=(14, height))
    ax_labels = figure.add_axes([0.05, 0.10, 0.18, 0.80])
    ax_pvals = figure.add_axes([0.22, 0.10, 0.12, 0.80])
    ax_ci = figure.add_axes([0.35, 0.10, 0.14, 0.80])
    ax_plot = figure.add_axes([0.51, 0.10, 0.44, 0.80])

    for axes in (ax_labels, ax_pvals, ax_ci):
        axes.axis("off")
        axes.set_ylim(-1, len(sequence))

    p_header = {"adj_pvalue": "Adjusted P-value",
                "p_value": "Nominal P-value"}.get(basis, "P-value")
    ax_labels.text(0, len(sequence), "Variant Call-type", weight="bold",
                   va="bottom", fontsize=13)
    ax_pvals.text(0.5, len(sequence), p_header, weight="bold",
                  va="bottom", ha="center", fontsize=13)
    ax_ci.text(0.5, len(sequence), "95% CI", weight="bold",
               va="bottom", ha="center", fontsize=13)

    for position, (kind, payload, colour) in enumerate(sequence):
        if kind == "spacer":
            continue
        if kind == "header":
            ax_labels.text(0, position, f"{payload} Genes", weight="bold",
                           color=colour, va="center", fontsize=12)
            ax_plot.axhline(position, color=colour, alpha=0.08, linewidth=12)
            continue

        row = payload
        # Bolding follows the basis the step chose.
        p_value = getattr(row, basis) if basis else float("nan")
        weight = "bold" if basis and p_value < ALPHA else "normal"

        ax_labels.text(0.1, position, row.a1type, va="center", fontsize=11,
                       fontweight=weight)
        shown = ("--" if basis is None
                 else (f"{p_value:.2e}" if p_value < 0.001 else f"{p_value:.3f}"))
        ax_pvals.text(0.5, position, shown, va="center", ha="center",
                      fontsize=11, fontweight=weight)
        ax_ci.text(0.5, position, f"[{row.ci_low:.2f}, {row.ci_high:.2f}]",
                   va="center", ha="center", fontsize=11, fontweight=weight)
        ax_plot.errorbar(
                row.odds_ratio, position,
                xerr=[[row.odds_ratio - row.ci_low],
                      [row.ci_high - row.odds_ratio]],
                fmt="o", color=colour, ecolor=colour, capsize=3, markersize=6,
                elinewidth=1.2, markeredgewidth=1.2)

    ax_plot.set_ylim(-1, len(sequence))
    ax_plot.axvline(1, color="black", linestyle="--", alpha=0.4, linewidth=1)
    ax_plot.set_xscale("log")
    ax_plot.set_xlabel(
            f"Odds ratio, unique sites in {case} relative to {reference} "
            f"(log scale, 95% CI)", fontsize=13, fontweight="bold")
    ax_plot.set_yticks([])

    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def main():
    args = parse_args()
    frame, case, reference, basis = load(args.input_dir)
    blocks = sections(frame, case, reference)
    if not blocks:
        raise SystemExit(
                f"no {SUBSET_SUFFIX.lstrip('_')} subset in {CSV_NAME} has a defined "
                f"odds ratio, so there is nothing to plot. Expected "
                f"'{case}{SUBSET_SUFFIX}' or '{reference}{SUBSET_SUFFIX}' in the "
                f"edit_dist column.")

    sequence = build_sequence(blocks)
    render(sequence, basis, case, reference, args.output, args.dpi)

    described = "nothing (no basis reached significance)" if basis is None else basis
    print(f"sections:   {', '.join(label for label, _, _ in blocks)}")
    print(f"call types: {sum(len(rows) for _, _, rows in blocks)}")
    print(f"bolded on:  {described}")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
