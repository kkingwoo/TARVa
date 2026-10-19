"""
dbpipe all: run the four database-build steps in their required order.

Each step remains the single writer of its own tables and runnable on its own.
The order is required: `load-calls` resolves sites to genes `annotate` wrote
and reads the cohort `load-samples` wrote, and `load-expression` keys on
`annotate`'s transcripts.

    `tarva dbpipe all --gtf <merged.gtf> --db <TARVa.db>
        --metadata <meta.csv> --rna-unique <rna_unique_dir>
        --stringtie-dir <assembly_dir> [every optional flag the four steps
        accept] [--overwrite]`

Input:  the GTF, the metadata CSV, the rna-unique VCF directory and the
        StringTie assembly directory. Arguments are the union of the four
        steps', so --db and --overwrite are given once; run
        `tarva dbpipe <step> --help` for which flags belong to which step.
Output: the database, written in this order:

            1. annotate         gene, transcript, exon
            2. load-samples     individual
            3. load-calls       site, call
            4. load-expression  expression
"""

from tarva import pipeline_runner
from tarva.dbpipe import annotate, load_calls, load_expression, load_samples

STEPS = (
    ("annotate", annotate),
    ("load-samples", load_samples),
    ("load-calls", load_calls),
    ("load-expression", load_expression),
    )


def add_arguments(parser):
    pipeline_runner.add_union(parser, STEPS)


def run(args):
    return pipeline_runner.run_steps(args, STEPS, "dbpipe")
