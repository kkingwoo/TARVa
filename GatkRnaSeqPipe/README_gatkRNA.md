README
================
kat j
2025-04-28

- [GATK RNAseq Pipeline](#gatk-rnaseq-pipeline)
  - [1. Sequencing read QC and
    alignment](#1-sequencing-read-qc-and-alignment)
  - [2. Alignment QC and
    variant-finding](#2-alignment-qc-and-variant-finding)

# GATK RNAseq Pipeline

This pipeline is last in the data pre-processing group of the workflow
and processes RNAseq data to create variant call format (vcf) files. The
following steps are executed for the pipeline:

## 1. Sequencing read QC and alignment

- [**Mono_MCI.slurm**](Mono_MCI.slurm) is the main slurm script for
  sequencing qc and alignment. Runs
  [**match_monocytes.py**](#match-monocytes).
  - [**match_monocytes.py**](match_monocytes.py)<a name="match-monocytes"></a>
    is the main python script for qc and alignment of sequencing data.
    1.  imports [**runqc.py**](runqc.py), which analyzes FastQC output
        and trims sequencing reads with Trimmomatic

    2.  creates and submits a slurm script to the job scheduler, which
        runs [**star.py**](#star), the main script that runs STAR in
        parallel for samples and imports [**star_run.py**](#star-run) as
        a custom module which executes STAR alignment for each sample.

## 2. Alignment QC and variant-finding

- [**full_gatk.slurm**](full_gatk.slurm) is the main slurm script for qc
  of alignments and variant-finding. Runs [**full_gatk.py**](#full-gatk)
  - [**full_gatk.py**](full_gatk.py)<a name="full-gatk"></a> runs the
    GATK pipeline in parallel and imports [**gatkPipe.py**]() as a
    custom module which contains commands for each step of the pipeline.
