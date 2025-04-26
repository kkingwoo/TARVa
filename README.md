README
================
kat j
2025-04-26

- [**T**ranscript **A**nalysis of **R**NA **Va**riants
  (TARVa)](#transcript-analysis-of-rna-variants-tarva)
  - [Project Overview](#project-overview)
  - [Tools, Environments and
    Dependencies](#tools-environments-and-dependencies)
  - [Directory Structure](#directory-structure)
  - [Workflow Overview](#workflow-overview)
  - [Limitations and Considerations](#limitations-and-considerations)
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
workflow is set up for use with paired-sample raw RNAseq and whole
genome sequencing files (fastq) from more than one condition and/or
tissue, and covers three data preprocessing pipelines in addition to the
main TARVa analysis pipeline.

## Tools, Environments and Dependencies

Python 3 is a requirement for this workflow and custom conda
environments were created for each group of pipelines ([data
preprocessing](#data-preprocessing-pipelines), [TARVa
analysis](#tarva-analysis-pipelines)).

For [**data preprocessing
pipelines**](#data-preprocessing-pipes)<a name="data-preprocessing-pipelines"></a>
([Stringtie](StringtiePipe/)<a name="stringtie-pipeline"></a>, [GATK
WGS](GatkWGSeqPipe/)<a name="gatk-wgs"></a>, [GATK
RNA](GatkRnaSeqPipe/))<a name="gatk-rna"></a>:  
- conda envs  
– *bio_qc* dependencies ([text](bio_qc_channels_dependencies.txt),
[yml](bio_qc.yml))

- R/4.3.3

- samtools/1.11

- gatk/4.2.0

- STAR

- Stringtie

- BWA

For the [**TARVa analysis
pipelines**](#tarva-analysis-pipes)<a name="tarva-analysis-pipelines"></a>
([Analysis set 1](TARVaCreation/OriginalBuild/), [Analysis set 1 checks
and changes](TARVaCreation/DownStreamDBCleanup/), [Analysis set
2](TARVaCreation/SecondSetBuild/)):  
- conda envs  
– *tarva* dependencies ([text](tarva_channels_dependencies.txt),
[yml](tarva.yml))

- python modules  
  – *sqlite3*  
  – *numpy*  
  – *scipy*  
  – *pandas*  
  – *statsmodels*  
  – *vcf*

## Directory Structure

<details>
<summary>
<strong>📂 Project Directory Tree (click to expand)</strong>
</summary>
<pre><code>.
├── GatkRnaSeqPipe
│   ├── README_gatkRNA.Rmd
│   ├── README_gatkRNA.md
│   ├── README_gatkRNA_files
│   │   └── figure-gfm
│   │       └── pressure-1.png
│   ├── STAR_array.slurm
│   ├── full_gatk.py
│   ├── full_gatk.slurm
│   ├── gatkPipe.py
│   ├── match_monocytes.py
│   ├── runqc.py
│   ├── star.py
│   └── star_run.py
├── GatkWGSeqPipe
│   ├── README_gatkWGS.Rmd
│   ├── README_gatkWGS.md
│   ├── README_gatkWGS_files
│   │   └── figure-gfm
│   │       └── pressure-1.png
│   ├── Step1_alignment.slurm
│   ├── Step2_mergeBams.slurm
│   ├── Step3_5_sortUnalnBams.slurm
│   ├── Step3_mkUnalnBam.slurm
│   ├── Step4_MergeBamAlns.py
│   ├── Step4_MergeBamAlns.slurm
│   ├── Step5_MD_SS.py
│   ├── Step5_MarkDuplicates_SortSam.slurm
│   ├── Step6.py
│   ├── Step6_BaseRecal_ApplyBQSR.slurm
│   ├── Step7.py
│   ├── Step7_CreatePons.py
│   ├── Step7_DBImport.py
│   ├── Step7_Mutect2_PON.slurm
│   ├── Step7_SelectVars.py
│   ├── Step8_ApplyPonsToAD.py
│   ├── Step8_Mutect2_AD.slurm
│   ├── Step9_GetPileupSummaries.slurm
│   ├── Step9_GetPileupSummaries_AD.py
│   └── Step9_GetPileupSummaries_Control.py
├── LICENSE
├── README.Rmd
├── README.md
├── StringtiePipe
│   ├── DGE_Pipe1_RunStar.slurm
│   ├── DGE_Pipe2_Round2.py
│   ├── DGE_Pipe2_StringRound1.slurm
│   ├── DGE_Pipe2_StringRound2.slurm
│   ├── DGE_Pipe2_StringRound3.slurm
│   ├── DGE_Pipe3_MakeBallgownFiles.py
│   ├── README_Stringtie.Rmd
│   ├── README_Stringtie.md
│   └── README_Stringtie_files
│       └── figure-gfm
│           └── pressure-1.png
├── TARVa
│   └── LICENSE
├── TARVaCreation
│   ├── DownStreamDBCleanup
│   │   ├── AnET.slurm
│   │   ├── README_set1_downstream.Rmd
│   │   ├── README_set1_downstream.md
│   │   ├── adjust_lists.py
│   │   ├── anet.py
│   │   ├── parse_vep.py
│   │   ├── sep_call_type.py
│   │   └── type_by_gene.py
│   ├── OriginalBuild
│   │   ├── BuildTarvaDBs.py
│   │   ├── README_OriginalBuild.Rmd
│   │   ├── README_OriginalBuild.md
│   │   ├── README_OriginalBuild_files
│   │   │   └── figure-gfm
│   │   │       └── pressure-1.png
│   │   ├── analyze_lens.py
│   │   ├── build.slurm
│   │   ├── checkRNA_againstWGS.py
│   │   ├── each_gene.py
│   │   ├── make_sample_tabs.py
│   │   ├── proc_pos.py
│   │   └── raw_counts.py
│   └── SecondSetBuild
│       ├── Build_SecondSet.py
│       ├── README_secondset_build.Rmd
│       ├── README_secondset_build.md
│       ├── Rel_Abu.slurm
│       ├── build2.slurm
│       ├── fpkm_within_condition_analysis.py
│       ├── get_second_set_counts.py
│       ├── make_secondset_tabs.py
│       ├── read_mstrg_tab.py
│       └── tpm_analysis.py
├── bio_qc.yml
├── bio_qc_channels_dependencies.txt
├── structure.txt
├── tarva.yml
└── tarva_channels_dependencies.txt
&#10;17 directories, 81 files
</code></pre>
</details>

## Workflow Overview

To implement the full TARVa workflow, the following [**data
pre-processing
pipelines**](#data-preprocessing-pipelines)<a name="data-preprocessing-pipes"></a>
are run in the following order:  
1. [Stringtie pipeline](StringtiePipe/)
([README](StringtiePipe/README_Stringtie.md))  
Stringtie receives STAR RNA sequencing alignment bam files as input and
creates transcript expression level files (.ctab) and a gtf file
containing nucleotide positions, transcript features, etc. of each
transcript identified in the dataset.

2.  [GATK WGS](GatkWGSeqPipe/)
    ([README](GatkWGSeqPipe/README_gatkWGS.md))  
    A GATK pipeline for small somatic variant-finding in WGS sequencing
    data is implemented to provide the TARVa tool with a baseline set of
    somatic genomic variants (vcf as output) for each sample. WGS data
    only for samples which also have RNAseq data available is used in
    this step.

3.  [GATK RNA](GatkRnaSeqPipe/)
    ([README](GatkRnaSeqPipe/README_gatkRNA.md))  
    For samples with WGS data that has been processed in the previous
    step, RNAseq sequencing data is run through a GATK pipeline for
    small variant finding in RNAseq (vcf as output).

Output from the data pre-processing phase of the workflow can then be
used as input for the [**TARVa analysis
pipelines**](#tarva-analysis-pipelines)<a name="tarva-analysis-pipes"></a>,
which can be executed in the following order:

1.  [Analysis set 1](TARVaCreation/OriginalBuild/)
    ([README](TARVaCreation/OriginalBuild/README_OriginalBuild.md))  
    RNA variants are first parsed to filter out any genomic variants,
    then a database is built from several data sources to store detailed
    information about RNA-specific modifications in the dataset.

2.  [Analysis set 1 checks and
    changes](TARVaCreation/DownStreamDBCleanup/)
    ([README](TARVaCreation/DownStreamDBCleanup/README_set1_downstream.md))  

3.  [Analysis set 2](TARVaCreation/SecondSetBuild/)
    ([README](TARVaCreation/SecondSetBuild/README_secondset_build.md))

## Limitations and Considerations

## License
