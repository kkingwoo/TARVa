"""
TARVa: command-line pipeline tool

Subparsers:

    - sort-input: sort raw FASTQs into the WGS|RNA/<group>/<adapter>/ tree the
        rest of the pipeline expects

        `tarva sort-input --metadata <meta.csv> --input <fastq_dir> --output <sorted_dir>`

    - prep-reference: decompress a reference, build its .fai/.dict, and build
        STAR and/or BWA indexes

        `tarva prep-reference --ref <genome.fa[.gz]> [--star-index DIR] [--bwa-index]`

    -qc: check read lengths, run fastQC and Trimmomatic, remove low quality tiles if needed
        ** Can run this more than once, until reads are satisfactorily cleaned for downstream analysis **

        `tarva qc --metadata <meta.csv> --input <fastq_dir> --output <output_dir>
        [--full] [--sample-size N] [--threads N] [--hpc <boolean>]
        [--partition <partition_name>] [--mail-user <mail-user@email>]`

    - stringtiepipe: first-pass STAR alignment, StringTie assembly in three
        rounds, and novel-transcript extraction. Runs before rnaseqpipe: the
        combined_reference.fa it writes is what the two-pass alignment maps to

        NOTE: `align` takes the .gtf; `assemble` and `novel-transcripts` need the
        .gff. StringTie and gffcompare share a parser that rejects the NCBI .gtf

        `tarva stringtiepipe align --metadata <meta.csv> --qc-output <qc_dir>
        --output <firstpass_dir> --ref <genome.fa> --star-index DIR
        --gtf <annotation.gtf> [--threads N] [--jobs N] [--hpc]`
        `tarva stringtiepipe assemble --metadata <meta.csv> --aligned <firstpass_dir>
        --output <assembly_dir> --gtf <annotation.gff>
        --strandedness rf|fr|unstranded [--stringtie-bin PATH] [--jobs N] [--hpc]`
        `tarva stringtiepipe novel-transcripts --assembly <assembly_dir>
        --gtf <annotation.gff> --ref <genome.fa> --output <novel_dir>
        [--gffcompare-bin PATH] [--gffread-bin PATH]`

    - wgseqpipe: WGS alignment, GATK preprocessing, panel of normals, variants

        `tarva wgseqpipe align --metadata <meta.csv> --qc-output <qc_dir>
        --output <aligned_dir> --ref <genome.fa> [--array] [--hpc]`
        `tarva wgseqpipe preprocess --metadata <meta.csv> --qc-output <qc_dir>
        --aligned <aligned_dir> --output <preprocessed_dir> --ref <genome.fa>
        [--known-sites <vcf>] [--array] [--hpc]`
        `tarva wgseqpipe pon --metadata <meta.csv> --preprocessed <preproc_dir>
        --output <pon_dir> --ref <genome.fa> --control-group <value>
        [--existing-pon <vcf>] [--array] [--hpc]`
        `tarva wgseqpipe variants --metadata <meta.csv> --preprocessed <preproc_dir>
        --output <vcf_dir> --ref <genome.fa> [--panel-of-normals <pon.vcf.gz>]
        [--array] [--hpc]`

    - rnaseqpipe: two-pass STAR alignment and GATK preprocessing for RNA
        variant calling; `preprocess` shares tarva/gatk_common.py with wgseqpipe

        `tarva rnaseqpipe align --metadata <meta.csv> --qc-output <qc_dir>
        --output <twopass_dir> --ref <combined_reference.fa> --star-index DIR [--array]`
        `tarva rnaseqpipe preprocess --metadata <meta.csv> --qc-output <qc_dir>
        --aligned <twopass_dir> --output <preprocessed_dir> --ref <combined_reference.fa>
        [--known-sites <vcf>] [--array] [--hpc]`
        `tarva rnaseqpipe variants --metadata <meta.csv> --preprocessed <preproc_dir>
        --output <vcf_dir> --ref <combined_reference.fa> [--array] [--hpc]`

    - analysispipe: pair the WGS and RNA call sets into the RNA-specific,
        panel-masked variant set the database build consumes

        `tarva analysispipe rna-unique --metadata <meta.csv>
        --wgs-variants <wgs_vcf_dir> --rna-variants <rna_vcf_dir>
        --output <rna_unique_dir> [--pon <pon.vcf.gz>]
        [--germline-resource <af-only.vcf.gz>] [--jobs N] [--overwrite]`
"""

import argparse
import sys

from tarva import sort_input as sort_input_module
from tarva import reference as reference_module
from tarva import qc as qc_module
from tarva.wgseqpipe import align as wgs_align
from tarva.wgseqpipe import preprocess as wgs_preprocess
from tarva.wgseqpipe import pon as wgs_pon
from tarva.wgseqpipe import variants as wgs_variants
from tarva.rnaseqpipe import align as rna_align
from tarva.rnaseqpipe import preprocess as rna_preprocess
from tarva.rnaseqpipe import variants as rna_variants
from tarva.stringtiepipe import align as stringtie_align
from tarva.stringtiepipe import assemble as stringtie_assemble
from tarva.stringtiepipe import novel_transcripts as stringtie_novel_transcripts
from tarva.analysispipe import rna_unique as analysis_rna_unique

WGSPIPE_STEPS = {
        "align": wgs_align,
        "preprocess": wgs_preprocess,
        "pon": wgs_pon,
        "variants": wgs_variants,
        }

RNASEQPIPE_STEPS = {
        "align": rna_align,
        "preprocess": rna_preprocess,
        "variants": rna_variants,
        }

STRINGTIEPIPE_STEPS = {
        "align": stringtie_align,
        "assemble": stringtie_assemble,
        "novel-transcripts": stringtie_novel_transcripts,
        }

ANALYSISPIPE_STEPS = {
        "rna-unique": analysis_rna_unique,
        }

def build_parser():
    parser = argparse.ArgumentParser(prog="tarva", description="TARVa multi-omics pipeline CLI")
    subparsers = parser.add_subparsers(dest="task", required=True)

    sort_input_parser = subparsers.add_parser(
            "sort-input", help=sort_input_module.__doc__.strip().splitlines()[0]
            )
    sort_input_module.add_arguments(sort_input_parser)
    sort_input_parser.set_defaults(func=sort_input_module.run)

    reference_parser = subparsers.add_parser(
            "prep-reference", help=reference_module.__doc__.strip().splitlines()[0]
            )
    reference_module.add_arguments(reference_parser)
    reference_parser.set_defaults(func=reference_module.run)

    qc_parser = subparsers.add_parser("qc", help=qc_module.__doc__.strip().splitlines()[0])
    qc_module.add_arguments(qc_parser)
    qc_parser.set_defaults(func=qc_module.run)

    stringtiepipe_parser = subparsers.add_parser(
            "stringtiepipe", help="StringTie assembly, including its first-pass STAR alignment"
            )
    stringtiepipe_sub = stringtiepipe_parser.add_subparsers(dest="step", required=True)
    for step_name, step_module in STRINGTIEPIPE_STEPS.items():
        step_parser = stringtiepipe_sub.add_parser(
                step_name, help=step_module.__doc__.strip().splitlines()[0]
                )
        step_module.add_arguments(step_parser)
        step_parser.set_defaults(func=step_module.run)

    wgseqpipe_parser = subparsers.add_parser("wgseqpipe", help="WGS data processing pipeline")
    wgseqpipe_sub = wgseqpipe_parser.add_subparsers(dest="step", required=True)
    for step_name, step_module in WGSPIPE_STEPS.items():
        step_parser = wgseqpipe_sub.add_parser(step_name, help=step_module.__doc__.strip().splitlines()[0])
        step_module.add_arguments(step_parser)
        step_parser.set_defaults(func=step_module.run)

    rnaseqpipe_parser = subparsers.add_parser(
            "rnaseqpipe", help="RNA-seq variant-calling pipeline"
            )
    rnaseqpipe_sub = rnaseqpipe_parser.add_subparsers(dest="step", required=True)
    for step_name, step_module in RNASEQPIPE_STEPS.items():
        step_parser = rnaseqpipe_sub.add_parser(
                step_name, help=step_module.__doc__.strip().splitlines()[0]
                )
        step_module.add_arguments(step_parser)
        step_parser.set_defaults(func=step_module.run)

    analysispipe_parser = subparsers.add_parser(
            "analysispipe", help="Pair WGS and RNA call sets into the RNA-specific variant set"
            )
    analysispipe_sub = analysispipe_parser.add_subparsers(dest="step", required=True)
    for step_name, step_module in ANALYSISPIPE_STEPS.items():
        step_parser = analysispipe_sub.add_parser(
                step_name, help=step_module.__doc__.strip().splitlines()[0]
                )
        step_module.add_arguments(step_parser)
        step_parser.set_defaults(func=step_module.run)

    return parser

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)

if __name__=="__main__":
    sys.exit(main())
