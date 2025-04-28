README
================
kat j
2025-04-28

- [Analysis set 2](#analysis-set-2)
  - [1. Addition of information to
    database](#1-addition-of-information-to-database)
  - [2. Relative Abundance](#2-relative-abundance)

# Analysis set 2

This is the last of the TARVa pipeline, adding variant call and sample
information for additional tissues and conditions. This step accepts as
input the RNA-seq variant information output (vcf) from the [**Gatk
RNAseq Pipeline**](../../GatkRnaSeqPipe/). Because modifications in the
additional samples will only be analyzed if they match a specified group
of modifications from the first analysis set, WGS data is not used.
Finally, assessment of relative abundance in samples is used to more
accurately inform the Variant Effect Prediction (VEP) analysis.

## 1. Addition of information to database

[**build2.slurm**](build2.slurm) submits
[**Build_SecondSet.py**](Build_SecondSet.py) for parsing and organizing
information from additional conditions and tissues, then populating the
database where there are matches with highlighted information from
Analysis Set 1. There are two @staticmethod scripts which
Build_SecondSet.py imports as custom modules:  
- [**make_secondset_tabs.py**](make_secondset_tabs.py) retrieves all
sample information and populates the database.  
- [**get_second_set_counts.py**](get_second_set_counts.py) calculates
counts for match instances of specified modifications in Analysis Set 1
also found in Analysis Set 2.

## 2. Relative Abundance

[**Rel_Abu.slurm**](Rel_Abu.slurm) submits
[**tpm_analysis.py**](tpm_analysis.py) and
[**fpkm_within_condition_analysis.py**](fpkm_within_condition_analysis.py)
which assess Stringtie TPM values for each sample and calculate FPKM
values per sample for genes and their transcripts within which specified
modifications occur. Both python scripts import read_mstrg_tab.py as a
module for parsing the Stringtie output data structure and files.
