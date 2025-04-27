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
be utilized for downstream analyses. [**Figure 1**](#figure) shows the
entity relationship diagram for the TARVa database.

![**Figure 1.**](../tarva_db_erd.png)<a name="figure"></a>  
**Figure 1.** ERD of TARVa database structure. *gtf_tab* contains
information from the merged Stringtie GTF, *fasta_tab* contains values
per nucleotide position in the reference fasta file, and *sample_tab*
contains values from individual VCF files for all conditions and
tissues. Yellow highlights: The reference nucleotide, alternate
nucleotide, and strand information were used to determine variant
call-type; Orange highlights: allele_depth, depth_reads, and length
retrieved to calculate the proportion of reads modified for each
instance of a variant call.

## 1. identify rna-specific variants

- [**checkRNA_againstWGS.py**](checkRNA_againstWGS.py) parses WGS and
  RNAseq vcf files for each sample, returning a new RNAseq vcf file
  containing only variants not found in the WGS data for the sample.

## 2. build database and populate tables

- [**build.slurm**]()

[**analyze_lens.py**](analyze_lens.py) creates a transcript length for
each gene, then calculates proportion of transcripts with a modification
for that gene. Due to the transient nature of RNA, the large number of
possible transcripts for some genes, the often-drastic length
differences among the transcripts, and the plausibility of differential
transcript usage and/or gene expression levels between conditions, a
normalization step was necessary. The type of information derived from
the current pipeline does not provide a clear-cut answer on which
specific transcript(s) the variant calls belong to. Thus, for each gene,
the minimum start and maximum end of all (start,end) positions for all
possible transcripts of that gene were used to assign a single
transcript length, per gene, for application in the normalization step.
The approach used here scales the information for all samples into the
space of n=1 transcript (and thus n=1 sequencing read) from which the
level of editing per gene, per sample, is derived.
