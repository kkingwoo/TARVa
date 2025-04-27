README
================
kat j
2025-04-27

- [GATK WGS Pipeline](#gatk-wgs-pipeline)
  - [Step 1: alignment](#step-1-alignment)
  - [Step 2: merge bams](#step-2-merge-bams)
  - [Step 3: make unaligned bams](#step-3-make-unaligned-bams)
  - [Step 3.5: sort unaligned bams](#step-35-sort-unaligned-bams)
  - [Step 4: merge bam alignments](#step-4-merge-bam-alignments)
  - [Step 5: mark duplicates](#step-5-mark-duplicates)
  - [Step 6: base recalibration and apply
    bqsr](#step-6-base-recalibration-and-apply-bqsr)
  - [Step 7: create panel of normals
    (pons)](#step-7-create-panel-of-normals-pons)
  - [Step 7: import genomics database](#step-7-import-genomics-database)
  - [Step 7: select variants](#step-7-select-variants)
  - [Step 8: variant calling in remaining or all
    conditions](#step-8-variant-calling-in-remaining-or-all-conditions)
  - [Step 9: filter vcf files](#step-9-filter-vcf-files)

# GATK WGS Pipeline

This pipeline is second in the data pre-processing group of the
workflow, and processes the whole genome sequencing (WGS) data to create
variant call format (vcf) files. The steps below are executed for the
pipeline, after FastQC and Trimmomatic have been applied for [**quality
checks and control of sequencing reads**](../GatkRnaSeqPipe/runqc.py).

## Step 1: alignment

- [**Step1_alignment.slurm**](Step1_alignment.slurm) runs BWA-mem to
  align the sequencing reads to a reference genome, then sorts the bams
  with samtools.

## Step 2: merge bams

- [**Step2_mergeBams.slurm**](Step2_mergeBams.slurm) uses samtools to
  merge sorted bams.

## Step 3: make unaligned bams

- [**Step3_mkUnalnBam.slurm**](Step3_mkUnalnBam.slurm) uses GATK’s
  *FastqToSam* to create unaligned bam files from the fastq files.

## Step 3.5: sort unaligned bams

- [**Step3_5_sortUnalnBams.slurm**](Step3_5_sortUnalnBams.slurm) uses
  GATK’s *SortSam* to sort the unaligned bams.

## Step 4: merge bam alignments

- [**Step4_MergeBamAlns.slurm**](Step4_MergeBamAlns.slurm) runs
  [**Step4_MergeBamAlns.py**](Step4_MergeBamAlns.py) which submits a
  slurm job that uses GATK’s *MergeBamAlignment* for merging aligned
  bams with unmapped bams.

## Step 5: mark duplicates

- [**Step5_MarkDuplicates_SortSam.slurm**](Step5_MarkDuplicates_SortSam.slurm)
  runs [**Step5_MD_SS.py**](Step5_MD_SS.py), which creates variable
  names then submits a job that runs GATK’s *MarkDuplicates* for marking
  duplicate alignments and *SortSam* to sort the marked duplicates
  output.

## Step 6: base recalibration and apply bqsr

- [**Step6_BaseRecal_ApplyBQSR.slurm**](Step6_BaseRecal_ApplyBQSR.slurm)
  submits [**Step6.py**](Step6.py), which submits GATK’s
  *BaseRecalibrator* and *ApplyBQSR* for improvement of accuracy of base
  quality estimates.

## Step 7: create panel of normals (pons)

This step is not needed if GATK’s downloadable PONs are sufficient for
the study.

- [**Step7_Mutect2_PON.slurm**](Step7_Mutect2_PON.slurm) is submitted to
  the job scheduler, which first runs [**Step7.py**](Step7.py) for
  Mutect2 variant calling of the Control samples in ‘tumor-only’ mode,
  then [**Step7_CreatePons.py**](Step7_CreatePons.py) which returns
  genomic coordinates and .

## Step 7: import genomics database

- [**Step7_DBImport.py**](Step7_DBImport.py) imports the pons to a
  database with GATK’s *GenomicsDBImport*.

## Step 7: select variants

- [**Step7_SelectVars.py**](Step7_SelectVars.py) runs GATK’s
  *SelectVars* with the previously created database as input, to create
  a final vcf file containing the set of PONs for downstream analysis.

## Step 8: variant calling in remaining or all conditions

- [**Step8_Mutect2_AD.slurm**](Step8_Mutect2_AD.slurm) runs
  [**Step8_ApplyPonsToAD.py**](Step8_ApplyPonsToAD.py) which uses a PONs
  as input to GATK’s *Mutect2* for variant-calling.

## Step 9: filter vcf files

- [**Step9_GetPileupSummaries.slurm**](Step9_GetPileupSummaries.slurm)
  runs separate python scripts for each condition
  ([**Step9_GetPileupSummaries_AD.py**](Step9_GetPileupSummaries_AD.py)
  and
  [**Step9_GetPileupSummaries_Control.py**](Step9_GetPileupSummaries_Control.py))
  for final filtration of the vcf files with GATK’s
  *GetPileupSummaries*, *CalculateContamination*, and
  *FilterMutectCalls*. The resulting files are used in the TARVa
  [**OriginalBuild**](../TARVaCreation/OriginalBuild/) pipeline.
