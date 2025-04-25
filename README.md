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
  - [Getting Started](#getting-started)
  - [How to Use](#how-to-use)
  - [Limitations and Considerations](#limitations-and-considerations)
  - [Acknowledgments](#acknowledgments)
  - [License](#license)

# **T**ranscript **A**nalysis of **R**NA **Va**riants (TARVa)

## Project Overview

[**TARVa**](#transcript-analysis-of-rna-variants-tarva) is a tool which
enables analysis of epitranscriptomic modifications. Starting with a
global overview, the pipeline moves through non-hierarchical
stratifications of the data, ending with local, per-position comparisons
across tissue types and conditions. The workflow is set up for use with
paired RNAseq and WGS raw sequencing datasets (fastq) from more than one
condition and/or tissue, and covers three data preprocessing pipelines
as well as the main TARVa pipeline.

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

For **data preprocessing pipelines** ([gatk_rna](GatkRnaSeqPipe/),
[gatk_wgs](GatkWGSeqPipe/), [stringtie](StringtiePipe/)):  
- conda envs - *bio_qc* ([text](bio_qc_channels_dependencies.txt),
[yml](bio_qc.yml))

## Getting Started

After installing the necessary tools and dependencies,

## How to Use

## Limitations and Considerations

## Acknowledgments

## License

\*\*\*\*\*\* NOTES NOTES
NOTES\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*
TARVa is a tool that identifies, analyzes, and compares RNA-editing
events between two conditions. This approach takes into consideration
both depth of read and length of transcript, and scans for “A-I” and
“C-U” variants in a strand-aware approach. TARVa also has an optional
feature, for studies which provide both RNA-seq and WGS-seq data for
each sample, to automatically filter out variant calls that are a result
of the DNA transcription process and not of post-modification
(RNA-editing) events.

###### **code draft suggestion 1; create user-input variable functions for paths, sample identifier information, and condition names**

###### **code draft suggestion 2: make necessary changes to reflect the optionality of genomic variant calls being added to the analysis pipeline**

#### Setup

See the “About.pdf” document for details on the . Not following this
pipeline might lead to inaccurate outputs from TARVa

##### 1.) Unzip TARVaScripts.zip

    Contents:  
        **a.** *BuildTarvaDBs.py*  
             This is the **main script**.
*make_sample_tabs*,*each_gene*, *raw_counts*, *dictionaries* and
*analyze_lens* are all static methods, housed also in the TARVaScripts
directory.  
        **b.** *TARVa.slurm*  
            This is the slurm directive script that activates the
**tarva** environment, assigns path information and runs the main
script.

**code draft suggestion 3: remove this from scripts directory and
instead provide information in a later section of this document
describing resources that are needed**

        **c.** *tarva_env_deps.txt*  
            TARVa runs in the custom conda environment, **tarva**. This
file indicates all of the packages and dependencies that belong in the
environment. To create the environment, [See step
2](#2-create-the-tarva-environment)

##### 2.) Create the Tarva environment

The **tarva** environment can be created using the command:  
`conda create -n [env_name] -c [channel] package(s)` c

##### 3.) Assign path names and condition names in scripts —–\>\>\>\>\> REMOVE THIS STEP ONCE code draft suggestion 1 is implemented

##### 4.) Submit job directive script
