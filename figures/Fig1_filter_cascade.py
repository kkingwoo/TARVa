"""
Figure 1: the filtering cascade, from loaded calls to the analysed gene set.

    python Fig1_filter_cascade.py --input-dir <dir> [--output Fig1.png]

Computes nothing: the counts and the drops come from the step. No main title.

Input:  filter_cascade.csv from `tarva analysispipe cascade`, in --input-dir.
Output: the image named by --output (the extension picks the format).
"""

import argparse

from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")

import matplotlib.pyplot as plt                    # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch   # noqa: E402

CSV_NAME = "filter_cascade.csv"

BOX_COLOUR = "#4878a8"
DROP_COLOUR = "#b03030"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True,
                        help="Directory holding filter_cascade.csv")
    parser.add_argument("--output", default="Fig1_filter_cascade.png",
                        help="Image to write (extension picks the format)")
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def load(input_dir):
    frame = pd.read_csv(Path(input_dir) / CSV_NAME)
    if frame.empty:
        raise SystemExit(f"{CSV_NAME} is empty -- run `tarva analysispipe cascade "
                         f"--output {input_dir}` first.")
    for column in ("stage", "description", "n_calls", "n_sites", "n_genes",
                   "n_calls_dropped"):
        if column not in frame.columns:
            raise SystemExit(f"{CSV_NAME} has no {column} column -- re-run cascade.")
    return frame


def render(frame, output, dpi):
    rows = list(frame.itertuples())
    widest = max(row.n_calls for row in rows) or 1

    figure, axes = plt.subplots(figsize=(9, 1.55 * len(rows) + 0.8))
    axes.set_xlim(0, 10)
    axes.set_ylim(0, 1.55 * len(rows))
    axes.axis("off")

    for index, row in enumerate(rows):
        # Top-down, so the first stage is the widest box at the top.
        top = 1.55 * (len(rows) - index)
        # Width in proportion to the call count, floored so a heavily
        # filtered final stage is still a readable box rather than a line.
        width = 2.2 + 6.6 * (row.n_calls / widest)
        left = 5.0 - width / 2

        axes.add_patch(FancyBboxPatch(
                (left, top - 1.0), width, 0.8,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                linewidth=1.2, edgecolor=BOX_COLOUR,
                facecolor=BOX_COLOUR, alpha=0.14))

        axes.text(5.0, top - 0.36, row.stage.replace("_", " "),
                  ha="center", va="center", fontsize=12, fontweight="bold")
        axes.text(5.0, top - 0.62,
                  f"{row.n_calls:,} calls   {row.n_sites:,} sites   "
                  f"{row.n_genes:,} genes",
                  ha="center", va="center", fontsize=10)
        axes.text(5.0 + width / 2 + 0.15, top - 0.5, row.description,
                  ha="left", va="center", fontsize=8.5, color="#555555",
                  style="italic")

        if index:
            axes.add_patch(FancyArrowPatch(
                    (5.0, top + 0.35), (5.0, top - 0.14),
                    arrowstyle="-|>", mutation_scale=14,
                    linewidth=1.1, color="#777777"))
            axes.text(5.25, top + 0.12, f"−{row.n_calls_dropped:,} calls",
                      ha="left", va="center", fontsize=9.5,
                      color=DROP_COLOUR, fontweight="bold")

    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def main():
    args = parse_args()
    frame = load(args.input_dir)
    render(frame, args.output, args.dpi)
    print(f"stages: {len(frame)}")
    print(f"from {int(frame['n_calls'].iloc[0]):,} calls to "
          f"{int(frame['n_calls'].iloc[-1]):,}")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
