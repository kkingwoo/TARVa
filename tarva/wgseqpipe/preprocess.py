"""
wgseqpipe preprocess: GATK preprocessing of aligned WGS BAMs

    `tarva wgseqpipe preprocess --metadata <meta.csv> --qc-output <qc_dir>
        --aligned <wgs_aligned_dir> --output <preprocessed_dir>
        --ref <genome.fa[.gz]> [--known-sites <vcf>] [--paired FWD,REV]
        [--temp-dir DIR] [--java-heap SIZE] [--keep-intermediates]
        [--read-group-repair auto|always|never] [--description TEXT]
        [--jobs N] [--sample-index N] [--overwrite] [--work-dir DIR]
        [--hpc] [--array] [--array-concurrency N] [--conda-env NAME]
        [--modules M1,M2] [--directives-dict "{...}"]`

Input is <aligned>/<group>/<individualID>/<individualID>_wgs.bam as written by
`wgseqpipe align`, plus the trimmed WGS FASTQ pair under --qc-output, which
FastqToSam needs for the unaligned BAM and the read-group header.

Output is <output>/<group>/<individualID>/<individualID>_bqsr.bam, or
<individualID>_coordsort.bam when --known-sites is omitted. Either way the BAM
is coordinate-sorted and indexed, ready for `wgseqpipe pon` and
`wgseqpipe variants`.
"""

from tarva.gatk_common import PipelineSpec, add_preprocess_arguments, run_preprocess
from tarva.wgseqpipe.align import ALIGNED_BAM_SUFFIX

SPEC = PipelineSpec(
        name="wgs",
        subcommand=("wgseqpipe", "preprocess"),
        is_rna=False,
        aligned_suffix=ALIGNED_BAM_SUFFIX,
        summary_name="wgs_preprocess_summary.csv",
        job_name="tarva_wgs_preprocess",
        )

REF_HELP = (
    "Reference genome FASTA, decompressed first if it ends in .gz. Must be the same "
    "reference `wgseqpipe align` used, and must already have a .fai and .dict beside "
    "it from `tarva prep-reference`"
    )

ALIGNED_HELP = (
    "The --output directory `wgseqpipe align` wrote; BAMs are read from "
    "<aligned>/<group>/<individualID>/<individualID>" + ALIGNED_BAM_SUFFIX
    )

KNOWN_SITES_HELP = (
    "Known-sites VCF for BaseRecalibrator. Optional: without it the chain stops after "
    "the coordinate sort and reports 'preprocessed (no BQSR)'. "
    )


def add_arguments(parser):
    add_preprocess_arguments(parser, SPEC, REF_HELP, ALIGNED_HELP, KNOWN_SITES_HELP)


def run(args):
    run_preprocess(args, SPEC)
