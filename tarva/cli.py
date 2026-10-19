"""
TARVa: a command-line workflow for analysing RNA variants in matched RNA-seq
and WGS data.

Four stages, run in order: read QC and data processing, the RNA-specific call
set, the database build, then the analyses. Every command and step carries its
own `--help`, which is where its flags and defaults are documented.

    `tarva <command> [<step>] [options]`
    `python -m tarva.cli <command> [<step>] [options]`

Input:  paired-end FASTQ files, a reference genome and annotation, and a
        metadata CSV naming each individual and its condition.
Output: one output directory per step, then a SQLite database and the
        analysis CSVs the figures read.

Commands:

    sort-input      sort raw FASTQs into the tree the rest of the workflow
                    expects
    prep-reference  decompress a reference, build its .fai/.dict, and build
                    STAR and/or BWA indexes
    qc              read lengths, FastQC and Trimmomatic; can be repeated
                    until the reads are clean enough for alignment
    stringtiepipe   align, assemble, novel-transcripts -- first-pass STAR
                    alignment, StringTie assembly in three rounds, and the
                    combined reference the two-pass alignment maps to
    wgseqpipe       align, preprocess, pon, variants -- WGS alignment, GATK
                    preprocessing, an optional panel of normals, and calls
    rnaseqpipe      align, preprocess, variants -- two-pass STAR alignment,
                    GATK preprocessing, and calls
    analysispipe    rna-unique for the RNA-specific call set, then the
                    analyses: cascade, gene-tests, threshold-sweep,
                    expression-correlation, site-uniqueness, local-sites, or
                    `all` for every analysis in order
    dbpipe          annotate, load-samples, load-calls, load-expression, or
                    `all` for the whole build in its required order

Three things that are easy to miss:

    - `stringtiepipe align` takes the .gtf, while `assemble` and
      `novel-transcripts` need the .gff -- StringTie and gffcompare share a
      parser that rejects an NCBI .gtf.
    - `stringtiepipe` runs before `rnaseqpipe`, which maps to the combined
      reference it writes.
    - `analysispipe rna-unique` produces what `dbpipe` consumes, so it is not
      part of `analysispipe all`.
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
from tarva.analysispipe import run_all as analysis_all
from tarva.analysispipe import cascade as analysis_cascade
from tarva.analysispipe import expression_correlation as analysis_expression
from tarva.analysispipe import gene_tests as analysis_gene_tests
from tarva.analysispipe import local_sites as analysis_local_sites
from tarva.analysispipe import site_uniqueness as analysis_site_uniqueness
from tarva.analysispipe import threshold_sweep as analysis_threshold_sweep
from tarva.dbpipe import run_all as db_all
from tarva.dbpipe import annotate as db_annotate
from tarva.dbpipe import load_samples as db_load_samples
from tarva.dbpipe import load_calls as db_load_calls
from tarva.dbpipe import load_expression as db_load_expression

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
        "all": analysis_all,
        "rna-unique": analysis_rna_unique,
        "cascade": analysis_cascade,
        "gene-tests": analysis_gene_tests,
        "threshold-sweep": analysis_threshold_sweep,
        "expression-correlation": analysis_expression,
        "site-uniqueness": analysis_site_uniqueness,
        "local-sites": analysis_local_sites,
        }

DBPIPE_STEPS = {
        "all": db_all,
        "annotate": db_annotate,
        "load-samples": db_load_samples,
        "load-calls": db_load_calls,
        "load-expression": db_load_expression,
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
            "analysispipe", help="The RNA-specific call set, then the analyses on the database"
            )
    analysispipe_sub = analysispipe_parser.add_subparsers(dest="step", required=True)
    for step_name, step_module in ANALYSISPIPE_STEPS.items():
        step_parser = analysispipe_sub.add_parser(
                step_name, help=step_module.__doc__.strip().splitlines()[0]
                )
        step_module.add_arguments(step_parser)
        step_parser.set_defaults(func=step_module.run)

    dbpipe_parser = subparsers.add_parser(
            "dbpipe", help="Build the TARVa database from the pipeline's output"
            )
    dbpipe_sub = dbpipe_parser.add_subparsers(dest="step", required=True)
    for step_name, step_module in DBPIPE_STEPS.items():
        step_parser = dbpipe_sub.add_parser(
                step_name, help=step_module.__doc__.strip().splitlines()[0]
                )
        step_module.add_arguments(step_parser)
        step_parser.set_defaults(func=step_module.run)

    return parser

def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    # Returned, not discarded: the `tarva` console script exits with this, so a
    # step's non-zero return reaches a shell or a scheduler instead of looking
    # like success. Steps that fail outright raise SystemExit with a message.
    return args.func(args)

if __name__=="__main__":
    sys.exit(main())
