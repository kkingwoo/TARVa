"""
analysispipe all: run every analysis step on a built database, in order.

`gene-tests` leads, because steps 2 to 5 read what it wrote -- `gene_class`
for the gene set and its subset labels, `run_meta` for the condition labels
and whether the set was filtered. `cascade` is independent of all of them but
runs last, since it reports the gene-level stages from `run_meta` rather than
recomputing them, and so prints a complete cascade only afterwards.

    `tarva analysispipe all --db <TARVa.db> [--output <dir>]
        [--min-depth 20] [--case-condition <label>]
        [--group-threshold 0.75] [--overwrite]`

Input:  a built database. Arguments are the union of the six steps',
        declared once each, so one --min-depth sets the eligibility rule
        everywhere and one --output collects every CSV. NOT included is
        `rna-unique`, which consumes VCF directories to build the call set
        the database is made FROM, and so runs before `dbpipe`.
Output: each step's own tables, CSVs and run_meta keys, in this order:

            1. gene-tests              the per-gene and per-call-type tests,
                                       the gene subsets, and each
                                       individual's call-type load
            2. threshold-sweep         sensitivity of the subsets to the
                                       enrichment threshold
            3. expression-correlation  expression against variant load
            4. site-uniqueness         unique versus common positions
            5. local-sites             the local level, per gene and
                                       call-type, then per position
            6. cascade                 the filtering cascade and background
                                       call-type rates
"""

from tarva import pipeline_runner
from tarva.analysispipe import (cascade, expression_correlation, gene_tests,
                                local_sites, site_uniqueness, threshold_sweep)

STEPS = (
    ("gene-tests", gene_tests),
    ("threshold-sweep", threshold_sweep),
    ("expression-correlation", expression_correlation),
    ("site-uniqueness", site_uniqueness),
    ("local-sites", local_sites),
    ("cascade", cascade),
    )


def add_arguments(parser):
    pipeline_runner.add_union(parser, STEPS)


def run(args):
    return pipeline_runner.run_steps(args, STEPS, "analysispipe")
