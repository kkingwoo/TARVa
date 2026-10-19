"""
analysispipe expression: gene-level expression matrices from StringTie ctab
    output, resolved through the step-1 locus map.

    `tarva analysispipe expression --stringtie <root> --db <tarva.db>
        --read-length 150 --output <results_dir>
        [--min-marker-fpkm 0.1] [--enst-only]
        [--hpc] [--r-module R/4.5.1] [--r-libs PATH] [--directives-dict "{...}"]`

    Input:  the StringTie ctab tree, the database for the cohort, and the maps
            `gtf-map` wrote into the same --output.
    Output: 02_db_samples.csv, 02_db_gene_map.csv, 02_sample_crosswalk.csv,
            02_gene_fpkm.csv, 02_gene_counts.csv, 02_gene_tpm.csv,
            02_gene_tmm.csv, 02_gene_tpm_all.csv, 02_gene_tmm_all.csv,
            02_excluded_fraction.csv, 02_marker_survival.csv,
            02_expression.rds, 02_expression_summary.json

    Python resolves the cohort and the directory-to-rids crosswalk from the
    database; R does the ctab ingest and matrix construction. The crosswalk is
    written to 02_sample_crosswalk.csv before R runs, so it can be checked.

    Genes are excluded by annotation, not by symbol: R reads the lists
    `gtf-map` wrote (00_excluded_genes.csv and 01_excluded_symbols.txt) and
    drops those genes before normalising, so they never enter the TPM
    denominator. Change what is excluded with `gtf-map --biotype-exclude` and
    re-run both steps -- the exclusion is not re-derived here, because it also
    determines which loci that step called chimeric.

    Two unfiltered matrices are written alongside the filtered ones for the
    confounding step, which models labelled genes filterByExpr had no reason
    to retain. They use the FILTERED per-sample denominators, so they are a
    strict superset of the filtered matrices -- identical values on shared
    genes -- rather than a second normalisation scale. TMM factors are
    estimated on the filtered matrix and carried over, since edgeR's trimmed
    mean is unreliable on a mostly near-zero matrix.
"""

import sys

from pathlib import Path

from tarva.analysispipe import db as tarva_db
from tarva.analysispipe.rstep import add_common_arguments, run_r_step


def add_arguments(parser):
    parser.add_argument(
            "--stringtie", required=True,
            help="Root directory searched recursively for t_data.ctab files",
            )
    parser.add_argument(
            "--db", required=True,
            help="TARVa SQLite database -- defines the cohort and supplies group labels",
            )
    parser.add_argument(
            "--read-length", type=int, required=True,
            help=(
                "Nominal read length for the prepDE-style count conversion "
                "(cov * length / read_length). Required rather than defaulted -- "
                "silently guessing it would be a bad failure mode. Affects only the "
                "count matrix feeding filterByExpr, not the TPM that feeds DSA."
                ),
            )
    parser.add_argument(
            "--table", default=tarva_db.DEFAULT_TABLE,
            help=f"Database table holding sample and gene info (default: {tarva_db.DEFAULT_TABLE})",
            )
    parser.add_argument(
            "--min-marker-fpkm", type=float, default=0.1,
            help=(
                "BRETIGEA markers with mean FPKM above this are kept regardless of "
                "filterByExpr (default: 0.1). Without this rescue, low-abundance "
                "microglial and endothelial markers are preferentially lost."
                ),
            )
    parser.add_argument(
            "--enst-only", action="store_true",
            help=(
                "Sensitivity option: keep only ENST transcripts, discarding novel "
                "isoforms within annotated loci instead of summing them into the parent gene"
                ),
            )
    add_common_arguments(parser)


def run(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not (output_dir / "01_locus_symbol_map.csv").is_file():
        raise SystemExit(
                f"tarva analysispipe: {output_dir / '01_locus_symbol_map.csv'} not found. "
                f"Run `tarva analysispipe gtf-map --output {output_dir}` first."
                )

    root = Path(args.stringtie)
    if not root.is_dir():
        raise SystemExit(f"tarva analysispipe: --stringtie is not a directory: {root}")

    samples = tarva_db.load_samples(args.db, args.table)
    samples.to_csv(output_dir / "02_db_samples.csv", index=False)
    counts = ", ".join(f"{g}={n}" for g, n in samples["con"].value_counts().items())
    print(f"Cohort from {args.table}: {len(samples)} samples ({counts})")

    tarva_db.load_gene_map(args.db, args.table).to_csv(
            output_dir / "02_db_gene_map.csv", index=False)

    crosswalk, missing = tarva_db.build_crosswalk(root, samples)
    if missing:
        for rid in missing:
            print(f"MISSING: no {tarva_db.DEFAULT_CTAB} at {root / rid}", file=sys.stderr)
        raise SystemExit(
                f"tarva analysispipe: {len(missing)} of {len(samples)} cohort sample(s) have no "
                f"{tarva_db.DEFAULT_CTAB} under {root}. Every sample in {args.table} must be "
                f"quantified, or composition estimates will not cover the cohort being analysed."
                )

    crosswalk_path = output_dir / "02_sample_crosswalk.csv"
    crosswalk.to_csv(crosswalk_path, index=False)
    print(f"Matched all {len(crosswalk)} cohort samples to a {tarva_db.DEFAULT_CTAB}")
    print(f"Wrote sample crosswalk: {crosswalk_path}")

    run_r_step(
            args,
            step="expression",
            script_name="02_build_matrix.R",
            extra={
                "stringtie_root": str(Path(args.stringtie).resolve()),
                "read_length": args.read_length,
                "min_marker_fpkm": args.min_marker_fpkm,
                # Nulled rather than dropped: params.json accumulates across
                # steps, so a value left by an older run would go on advertising
                # a symbol-regex filter that no longer exists.
                "exclude_pattern": None,
                "enst_only": bool(args.enst_only),
                "db_table": args.table,
                },
            )
