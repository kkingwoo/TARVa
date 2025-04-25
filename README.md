---
title: "README"
author: "kat j"
date: "2025-04-25"
output: 
  github_document:
    toc: true
    toc_depth: 3
    number_sections: false
    preserve_yaml: true
---

README
================
kat j
2025-04-25

- [**T**ranscript **A**nalysis of **R**NA **Va**riants
  (TARVa)](#transcript-analysis-of-rna-variants-tarva)
  - [Project Overview](#project-overview)
  - [Getting Started](#getting-started)
  - [Environments and Dependencies](#environments-and-dependencies)
  - [How to Use](#how-to-use)
  - [Limitations and Considerations](#limitations-and-considerations)
  - [Acknowledgments](#acknowledgments)
  - [License](#license)

# **T**ranscript **A**nalysis of **R**NA **Va**riants (TARVa)

## Project Overview

## Getting Started

## Environments and Dependencies

## How to Use

## Limitations and Considerations

## Acknowledgments

## License

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
