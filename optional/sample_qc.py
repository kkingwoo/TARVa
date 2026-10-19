"""
analysispipe sample-qc: sex check, PCA and correlation outlier detection on the
    expression matrix, before anything is deconvolved.

    `tarva analysispipe sample-qc --db <tarva.db> --output <results_dir>
        [--sex-genes XIST,RPS4Y1,DDX3Y,UTY,KDM5D] [--n-variable 2000]
        [--cor-sd 3] [--covariate-file <pmi.csv>] [--covariate-file-key individualID]
        [--hpc] [--r-module R/4.5.1] [--r-libs PATH] [--directives-dict "{...}"]`

    Input:  the matrices `expression` wrote into the same --output, and the
            database for the recorded sex and other covariates.
    Output: 03_covariates.csv, 03_sample_qc.csv, 03_pc_associations.csv,
            03_sample_qc_summary.json, 03_pca.png

    Nothing is excluded automatically. Proposed exclusions are reported and
    written to 03_sample_qc.csv for a decision -- dropping a sample changes
    the cohort, so the call belongs to the user rather than the pipeline.

    The sex check does not assume how msex is coded. Samples are split into two
    groups by their X/Y expression contrast and cross-tabulated against the
    recorded value; a clean 2x2 establishes the coding, and mismatches are the
    off-diagonal cells.

    PC1 to PC3 are correlated against the non-mRNA fraction step 2 removed and
    against library size. That table is how the annotation-based exclusion is
    judged: if the leading component still tracks how much was taken out, the
    matrix is still dominated by library composition rather than by biology.
"""

from pathlib import Path

from tarva.analysispipe import db as tarva_db
from tarva.analysispipe.rstep import add_common_arguments, run_r_step

DEFAULT_SEX_GENES = "XIST,RPS4Y1,DDX3Y,UTY,KDM5D"


def _csv_list(raw):
    return [v.strip() for v in str(raw).split(",") if v.strip()]


def add_arguments(parser):
    parser.add_argument(
            "--db", required=True,
            help="TARVa SQLite database -- supplies msex and the other covariates",
            )
    parser.add_argument(
            "--table", default=tarva_db.DEFAULT_TABLE,
            help=f"Database table holding sample info (default: {tarva_db.DEFAULT_TABLE})",
            )
    parser.add_argument(
            "--sex-genes", default=DEFAULT_SEX_GENES,
            help=(
                "Comma-separated sex-informative genes (default: "
                f"{DEFAULT_SEX_GENES}). XIST is treated as X-linked, the rest as Y-linked."
                ),
            )
    parser.add_argument(
            "--n-variable", type=int, default=2000,
            help="Most variable genes used for PCA and sample correlation (default: 2000)",
            )
    parser.add_argument(
            "--cor-sd", type=float, default=3.0,
            help=(
                "Flag a sample whose median correlation to the others is this many "
                "MAD-scaled deviations below the cohort median (default: 3)"
                ),
            )
    tarva_db.add_covariate_file_arguments(parser)
    add_common_arguments(parser)


def run(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not (output_dir / "02_expression.rds").is_file():
        raise SystemExit(
                f"tarva analysispipe: {output_dir / '02_expression.rds'} not found. "
                f"Run `tarva analysispipe expression --output {output_dir}` first."
                )

    covariates = tarva_db.load_covariates(args.db, args.table)
    if args.covariate_file:
        covariates = tarva_db.merge_covariate_file(
                covariates, args.covariate_file,
                key=args.covariate_file_key,
                columns=_csv_list(args.covariate_file_columns) if args.covariate_file_columns else None,
                )
    covariates.to_csv(output_dir / "03_covariates.csv", index=False)
    print(f"Loaded covariates for {len(covariates)} samples: "
          + ", ".join(c for c in covariates.columns if c not in ("rids", "con")))

    run_r_step(
            args,
            step="sample_qc",
            script_name="03_sample_qc.R",
            extra={
                "sex_genes": _csv_list(args.sex_genes),
                "n_variable": args.n_variable,
                "cor_sd": args.cor_sd,
                },
            )
