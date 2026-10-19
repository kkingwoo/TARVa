"""
analysispipe gtf-map: build the locus->symbol map from the merged StringTie GTF
    and decide whether that merged assembly is usable for deconvolution.

    `tarva analysispipe gtf-map --gtf <merged.gtf> --output <results_dir>
        [--biotype-tsv <annotation.tsv>] [--biotype-exclude Mt,rRNA,tRNA]
        [--biotype-id-column geneID] [--biotype-type-column transcriptType]
        [--cell-types ast,end,...] [--gate-cell-types ast,end,mic,neu,oli]
        [--top-n 500] [--gate-clean-min 300] [--gate-clean-floor 200]
        [--hpc] [--r-module R/4.5.1] [--r-libs PATH] [--directives-dict "{...}"]`

    This is a blocking gate -- run it before anything else.

    Input:  the merged StringTie GTF, and optionally an annotation TSV of gene
            biotypes for the exclusions below.
    Output: into --output -- 00_excluded_genes.csv, 01_locus_symbol_map.csv,
            01_chimeric_symbols.txt, 01_excluded_symbols.txt,
            01_fragmented_symbols.csv, 01_marker_coverage.csv,
            01_gtf_map_summary.json

    It is also where genes are excluded by annotation, because the chimerism
    test below is what they distort: a mitochondrial or rRNA gene overlapping a
    protein-coding one makes the merged assembly fuse their loci, and the gate
    then blacklists the protein-coding gene for a collision with something that
    was never going to be analysed. Removing them first is not the same filter
    applied earlier -- it changes which loci are called chimeric at all.

    Quantification used a merged assembly, so t_data.ctab's gene_name column is
    "." for novel isoforms inside annotated loci. Joining ctab rows to symbols
    through gene_id (the MSTRG locus) instead keeps that expression, which is
    what this map is for.

    Two failure modes are measured:
        fragmentation - one symbol spread across several MSTRG loci. Harmless;
                        step 2 aggregates on symbol and puts it back together.
        chimerism     - one MSTRG locus absorbing two or more reference genes.
                        Not fixable, since expression at that locus cannot be
                        attributed to either gene. Those symbols are blacklisted,
                        which is affordable because BRETIGEA supplies 1000
                        candidate markers per cell type and only 100 are used.
"""

from pathlib import Path

from tarva.analysispipe import biotypes
from tarva.analysispipe.rstep import add_common_arguments, run_r_step

DEFAULT_GATE_CELL_TYPES = "ast,end,mic,neu,oli"


def _csv_list(raw):
    return [v.strip() for v in str(raw).split(",") if v.strip()]


def add_arguments(parser):
    parser.add_argument(
            "--gtf", required=True,
            help="Merged StringTie GTF (the one passed to stringtie -e -B -G). May be gzipped.",
            )
    parser.add_argument(
            "--cell-types", default=None,
            help=(
                "Comma-separated cell types to report marker coverage for. Default: every "
                "cell type present in BRETIGEA::markers_df_human_brain, so OPC is included "
                "here for the record even though it is excluded downstream."
                ),
            )
    parser.add_argument(
            "--gate-cell-types", default=DEFAULT_GATE_CELL_TYPES,
            help=(
                "Comma-separated cell types the pass/fail gate applies to "
                f"(default: {DEFAULT_GATE_CELL_TYPES}) -- the five carried forward "
                "following Wang et al. 2020, which excluded OPC for marker incoherence."
                ),
            )
    parser.add_argument(
            "--top-n", type=int, default=500,
            help="Rank-ordered markers per cell type to consider, per Pottier et al. (default: 500)",
            )
    parser.add_argument(
            "--gate-clean-min", type=int, default=300,
            help="Every gated cell type needs at least this many clean markers to pass (default: 300)",
            )
    parser.add_argument(
            "--gate-clean-floor", type=int, default=200,
            help=(
                "Below this many clean markers for any gated cell type the merged assembly "
                "is unusable and the step exits non-zero (default: 200)"
                ),
            )
    biotypes.add_biotype_arguments(parser)
    add_common_arguments(parser)


def run(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    exclusions = biotypes.resolve_exclusions(args, output_dir)

    run_r_step(
            args,
            step="gtf_map",
            script_name="01_gtf_map.R",
            extra={
                "merged_gtf": args.gtf,
                "cell_types": _csv_list(args.cell_types) if args.cell_types else None,
                "gate_cell_types": _csv_list(args.gate_cell_types),
                "top_n": args.top_n,
                "gate_clean_min": args.gate_clean_min,
                "gate_clean_floor": args.gate_clean_floor,
                "biotype_tsv": args.biotype_tsv,
                "biotype_exclude": biotypes.csv_list(args.biotype_exclude),
                "biotype_exclusions": str(exclusions) if exclusions else None,
                },
            )
