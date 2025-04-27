README
================
kat j
2025-04-27

- [Building the TARVa database](#building-the-tarva-database)
  - [1. identify rna-specific
    variants](#1-identify-rna-specific-variants)
  - [2. build database and populate
    tables](#2-build-database-and-populate-tables)

# Building the TARVa database

After data pre-processing pipelines are complete, information from vcfs,
fasta, gtf, and other sources are integrated into a database which can
be utilized for downstream analyses.

## 1. identify rna-specific variants

- [**checkRNA_againstWGS.py**](checkRNA_againstWGS.py) parses WGS and
  RNAseq vcf files for each sample, returning a new RNAseq vcf file
  containing only variants not found in the WGS data for the sample.

## 2. build database and populate tables

- [**build.slurm**](build.slurm) submits
  [**BuildTarvaDBs.py**](BuildTarvaDBs.py) for performing various
  operations on the data and building and populating the database. The
  following scripts are @staticmethods which are imported into
  BuildTarvaDBs.py as modules:
  - [**dictionaries.py**](dictionaries.py) has multiple functions
    related to populating dictionaries

  - [**make_sample_tabs.py**](make_sample_tabs.py) which organizes data
    related to each sample, for populating the **sample_tab** table in
    the database.  

  - [**analyze_lens.py**](analyze_lens.py) creates a transcript length
    for each gene and then for each modification, calculates proportion
    of transcripts for the gene that have that modification. Due to the
    transient nature of RNA, the large number of possible transcripts
    for some genes, the often-drastic length differences among the
    transcripts, and the plausibility of differential transcript usage
    and/or gene expression levels between conditions, a normalization
    step was necessary. The type of information derived from the current
    pipeline does not provide a clear-cut answer on which specific
    transcript(s) the variant calls belong to. Thus, for each gene, the
    minimum start and maximum end of all (start,end) positions for all
    possible transcripts of that gene were used to assign a single
    transcript length, per gene, for application in the normalization
    step. The approach used here scales the information for all samples
    into the space of n=1 transcript (and thus n=1 sequencing read) from
    which the level of editing per gene, per sample, is derived. The
    equation and variables are described below:
