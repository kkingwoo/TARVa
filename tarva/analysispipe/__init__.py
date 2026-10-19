"""
analysispipe: the RNA-specific call set, the analyses that read the database,
and the optional cell-type deconvolution.

Three bodies of work live here.

**rna-unique** pairs each sample's WGS and RNA VCFs and writes the
RNA-specific, panel-masked variant set the database build consumes.

**The analyses** -- cascade, gene_tests, threshold_sweep,
expression_correlation, site_uniqueness and local_sites, with `run_all` to run
them in order -- read the built database and write their tables and CSVs.
significance.py holds the one p-value convention they share.

**The cell-type deconvolution** -- gtf_map, expression, sample_qc, coherence,
dsa, stats, validate, confounding -- plus external_checks, which are OPTIONAL
and are not currently registered in cli.py. Registering one means re-adding
its import and its ANALYSISPIPE_STEPS entry there; nothing in the modules has
to change. Each step is separately runnable, writes into a shared --output
directory so intermediate results can be inspected, and is numbered in the
order it must run. Those analyses are done in R -- see the R/ subdirectory --
with Python handling parameter resolution, SLURM submission and all database
access; rstep.py wires the two halves together.

The deconvolution follows the Digital Sorting Algorithm of Zhong et al. 2013
with the BRETIGEA marker sets of McKenzie et al. 2018, as applied in Wang et
al. 2020 (Mol Neurodegener 15:38) and Pottier et al. 2022 (Brain 145:2472).
"""
