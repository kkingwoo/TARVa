README
================
kat j
2025-04-25

- [**T**ranscript **A**nalysis of **R**NA **Va**riants
  (TARVa)](#transcript-analysis-of-rna-variants-tarva)
  - [Project Overview](#project-overview)
  - [Directory Structure](#directory-structure)
  - [Tools, Environments and
    Dependencies](#tools-environments-and-dependencies)
  - [Workflow Overview](#workflow-overview)
  - [License](#license)

# **T**ranscript **A**nalysis of **R**NA **Va**riants (TARVa)

## Project Overview

[**TARVa**](#transcript-analysis-of-rna-variants-tarva) is a tool which
enables analysis of epitranscriptomic modifications from raw RNA
sequencing data, taking into consideration both depth of read and length
of transcript, and scans for all variant call-types in a strand-aware
fashion. Starting with a global overview, the pipeline moves through
non-hierarchical stratifications of the data, ending with local,
per-position comparisons across tissue types and conditions. The
workflow is set up for use with paired RNAseq and WGS raw sequencing
datasets (fastq) from more than one condition and/or tissue, and covers
three data preprocessing pipelines as well as the main TARVa pipeline.

## Directory Structure

<details>
<summary>
<strong>📂 Project Directory Tree (click to expand)</strong>
</summary>
<pre><code>.
├── TARVa
│   ├── GatkRnaSeqPipe
│   │   ├── README_gatkRNA.Rmd
│   │   ├── full_gatk.py
│   │   ├── full_gatk.slurm
│   │   └── gatkPipe.py
│   ├── GatkWGSeqPipe
│   │   ├── README_gatkWGS.Rmd
│   │   ├── Step1_alignment.slurm
│   │   ├── Step2_mergeBams.slurm
│   │   ├── Step3_5_sortUnalnBams.slurm
│   │   ├── Step3_mkUnalnBam.slurm
│   │   ├── Step4_MergeBamAlns.py
│   │   ├── Step4_MergeBamAlns.slurm
│   │   ├── Step5_MD_SS.py
│   │   ├── Step5_MarkDuplicates_SortSam.slurm
│   │   ├── Step6.py
│   │   ├── Step6_BaseRecal_ApplyBQSR.slurm
│   │   ├── Step7.py
│   │   ├── Step7_CreatePons.py
│   │   ├── Step7_DBImport.py
│   │   ├── Step7_Mutect2_PON.slurm
│   │   ├── Step7_SelectVars.py
│   │   ├── Step8_ApplyPonsToAD.py
│   │   ├── Step8_Mutect2_AD.slurm
│   │   ├── Step9_GetPileupSummaries.slurm
│   │   ├── Step9_GetPileupSummaries_AD.py
│   │   └── Step9_GetPileupSummaries_Control.py
│   ├── LICENSE
│   ├── README.Rmd
│   ├── README.md
│   ├── StringtiePipe
│   │   ├── DGE_Pipe1_RunStar.slurm
│   │   ├── DGE_Pipe2_Round2.py
│   │   ├── DGE_Pipe2_StringRound1.slurm
│   │   ├── DGE_Pipe2_StringRound2.slurm
│   │   ├── DGE_Pipe2_StringRound3.slurm
│   │   ├── DGE_Pipe3_MakeBallgownFiles.py
│   │   └── README_Stringtie.Rmd
│   ├── TARVa
│   │   └── LICENSE
│   ├── TARVaCreation
│   │   ├── DownStreamDBCleanup
│   │   │   ├── AnET.slurm
│   │   │   ├── adjust_lists.py
│   │   │   ├── anet.py
│   │   │   ├── parse_vep.py
│   │   │   ├── sep_call_type.py
│   │   │   └── type_by_gene.py
│   │   ├── OriginalBuild
│   │   │   ├── BuildTarvaDBs.py
│   │   │   ├── analyze_lens.py
│   │   │   ├── build.slurm
│   │   │   ├── checkRNA_againstWGS.py
│   │   │   ├── each_gene.py
│   │   │   ├── make_sample_tabs.py
│   │   │   ├── proc_pos.py
│   │   │   └── raw_counts.py
│   │   └── SecondSetBuild
│   │       ├── Build_SecondSet.py
│   │       ├── Rel_Abu.slurm
│   │       ├── build2.slurm
│   │       ├── fpkm_within_condition_analysis.py
│   │       ├── get_second_set_counts.py
│   │       ├── make_secondset_tabs.py
│   │       ├── read_mstrg_tab.py
│   │       └── tpm_analysis.py
│   ├── bio_qc.yml
│   ├── bio_qc_channels_dependencies.txt
│   ├── tarva.yml
│   └── tarva_channels_dependencies.txt
└── structure.txt
&#10;10 directories, 63 files
</code></pre>
</details>

## Tools, Environments and Dependencies

Python 3 is a requirement for this workflow.

For **data preprocessing pipelines** ([gatk_rna](GatkRnaSeqPipe/),
[gatk_wgs](GatkWGSeqPipe/), [stringtie](StringtiePipe/)):  
- conda envs  
– *bio_qc* ([text](bio_qc_channels_dependencies.txt), [yml](bio_qc.yml))

- R/4.3.3

- samtools/1.11

- gatk/4.2.0

- STAR

- Stringtie

- BWA

For the **TARVa analysis pipelines** ([analysis set
1](TARVaCreation/OriginalBuild/), [analysis set 1: checks and
changes](TARVaCreation/DownStreamDBCleanup/), [analysis set
2](TARVaCreation/SecondSetBuild/)):  
- conda envs  
– *tarva* ([text](tarva_channels_dependencies.txt), [yml](tarva.yml))

- python modules  
  – *sqlite3*  
  – *numpy*  
  – *scipy*  
  – *pandas*  
  – *statsmodels*  
  – *vcf*

## Workflow Overview

To implement the full TARVa workflow, the **data pre-processing
pipelines** ( [stringtie](StringtiePipe/)
([README](StringtiePipe/README_Stringtie.md)),
[gatk_rna](GatkRnaSeqPipe/)
([README](GatkRnaSeqPipe/README_gatkRNA.md)), [gatk_wgs](GatkWGSeqPipe/)
([README](GatkWGSeqPipe/README_gatkWGS.md))) are run first.  
Output from each pipeline can then be used as input for the **TARVa
analysis pipeline** ([analysis set 1](TARVaCreation/OriginalBuild/)
([README](TARVaCreation/OriginalBuild/README_OriginalBuild.md)),
[analysis set 1: checks and changes](TARVaCreation/DownStreamDBCleanup/)
([README](TARVaCreation/DownStreamDBCleanup/README_set1_downstream.md)),
[analysis set 2](TARVaCreation/SecondSetBuild/)
([README](TARVaCreation/SecondSetBuild/README_secondset_build.md))).

## License
