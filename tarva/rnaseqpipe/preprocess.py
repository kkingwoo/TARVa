"""
rnaseqpipe preprocess: GATK preprocessing of two-pass-aligned RNA BAMs

    `tarva rnaseqpipe preprocess --metadata <meta.csv> --qc-output <qc_dir>
        --aligned <rna_aligned_dir> --output <preprocessed_dir>
        --ref <genome.fa[.gz]> [--known-sites <vcf>] [--paired FWD,REV]
        [--temp-dir DIR] [--java-heap SIZE] [--keep-intermediates]
        [--read-group-repair auto|always|never] [--description TEXT]
        [--jobs N] [--sample-index N] [--overwrite] [--work-dir DIR]
        [--hpc] [--array] [--array-concurrency N] [--conda-env NAME]
        [--modules M1,M2] [--directives-dict "{...}"]`

Input is the SECOND STAR alignment -- `rnaseqpipe align`'s `--twopassMode Basic`
run for variant calling -- not the first-pass BAM that StringTie used for alignment. Those
are different alignments of the same reads for different purposes. 
See documentation for more details.

"""

from tarva.gatk_common import PipelineSpec, add_preprocess_arguments, run_preprocess
from tarva.rnaseqpipe.align import TWOPASS_BAM_SUFFIX

SPEC = PipelineSpec(
        name="rna",
        subcommand=("rnaseqpipe", "preprocess"),
        is_rna=True,
        aligned_suffix=TWOPASS_BAM_SUFFIX,
        summary_name="rna_preprocess_summary.csv",
        job_name="tarva_rna_preprocess",
        )

REF_HELP = (
    "Reference genome FASTA, decompressed first if it ends in .gz. Normally the "
    "combined_reference.fa from `stringtiepipe novel-transcripts`, not the plain "
    "genome -- it must be the same reference the two-pass STAR alignment used, and must "
    "already have a .fai and .dict beside it from `tarva prep-reference`"
    )

ALIGNED_HELP = (
    "The --output directory `rnaseqpipe align` wrote -- the two-pass STAR run for "
    "variant calling, NOT `stringtiepipe align`'s first pass. BAMs are read from "
    "<aligned>/<group>/<individualID>/<individualID>" + TWOPASS_BAM_SUFFIX
    )

KNOWN_SITES_HELP = (
    "Known-sites VCF for BaseRecalibrator. Optional: without it the chain stops after "
    "SplitNCigarReads (and AddOrReplaceReadGroups) and reports 'preprocessed (no BQSR)'. "
    "Contig names must match --ref, which for the combined reference means the novel "
    "transcript contigs are simply absent from the VCF -- BQSR ignores "
    "contigs it has no sites for"
    )


def add_arguments(parser):
    add_preprocess_arguments(parser, SPEC, REF_HELP, ALIGNED_HELP, KNOWN_SITES_HELP)


def run(args):
    run_preprocess(args, SPEC)
