"""
TARVa: command-line pipeline tool

Subparsers:
    -qc: check read lengths, run fastQC and Trimmomatic, remove low quality tiles if needed
        ** Can run this more than once, until reads are satisfactorily cleaned for downstream analysis **
        
        `tarva qc --metadata <meta.csv> --input <fastq_dir> --output <output_dir> 
        [--full] [--sample-size N] [--threads N] [--hpc <boolean>]
        [--partition <partition_name>] [--mail-user <mail-user@email>]`

    - wgseqpipe: data processing and QC for raw WGS data

        `tarva wgspipe align --metadata <meta.csv> --wgs <wgs_fastq_path> --output <output_path> 
        --ref <genomic_ref.fa> --paired <fw_pattern>,<rev_pattern> --hpc <boolean> 
        --partition <partition_name> --mail-user <mail-user@email>`

    - rnaseqpipe: data processing and QC for raw RNA-seq data

    - stringtiepipe: run StringTie on RNA-seq data

    - analysispipe: build local database and analyze processed reads
"""

import argparse
import sys

from tarva import qc as qc_module
from tarva.wgseqpipe import align as wgs_align

WGSPIPE_STEPS = {
        "align": wgs_align,
        }

def build_parser():
    parser = argparse.ArgumentParser(prog="tarva", description="TARVa multi-omics pipeline CLI")
    subparsers = parser.add_subparsers(dest="task", required=True)

    qc_parser = subparsers.add_parser("qc", help=qc_module.__doc__.strip().splitlines()[0])
    qc_module.add_arguments(qc_parser)
    qc_parser.set_defaults(func=qc_module.run)

    wgseqpipe_parser = subparsers.add_parser("wgseqpipe", help="WGS data processing pipeline")
    wgseqpipe_sub = wgseqpipe_parser.add_subparsers(dest="step", required=True)
    for step_name, step_module in WGSPIPE_STEPS.items():
        step_parser = wgseqpipe_sub.add_parser(step_name, help=step_module.__doc__.strip().splitlines()[0])
        step_module.add_arguments(step_parser)
        step_parser.set_defaults(func=step_module.run)

    # rnaseqpipe

    # stringtiepipe

    # analysispipe

    return parser

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)

if __name__=="__main__":
    sys.exit(main())


