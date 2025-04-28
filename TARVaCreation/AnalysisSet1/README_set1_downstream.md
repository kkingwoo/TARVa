README
================
kat j
2025-04-28

- [Analysis set 1](#analysis-set-1)
  - [1. global level distribution assignment of modified
    positions](#1-global-level-distribution-assignment-of-modified-positions)
  - [2. conduct global and local-level
    analyses](#2-conduct-global-and-local-level-analyses)
  - [3. Plot data](#3-plot-data)

# Analysis set 1

This analysis step is part of the TARVa analysis pipeline, and is run
after the database has been built. Operations such as accessing
additional information about the genes, needed for conducting the global
and local analyses are carried out in this section.

## 1. global level distribution assignment of modified positions

- [**AnET.slurm**](AnET.slurm) submits [**anet.py**](anet.py) to the job
  scheduler to carry out operations such as retrieving additional gene
  information and global-level distribution assignment of genes.

  - [**type_by_gene.py**](type_by_gene.py) is imported to anet.py as a
    custom module. This module carries out assessment and assignment of
    variant call-types in a strand-aware fashion, and analysis and
    assignment of modified positions to their respective global-level
    distribution subsets for each gene: (*common*=commonly modified and
    *unique*=uniquely modified)

## 2. conduct global and local-level analyses

- [**Downstream.slurm**](Downstream.slurm) submits
  [**downstream_analysis.py**](downstream_analysis.py) to carry out the
  remainder of global and local-level analyses.

## 3. Plot data

[**Figures.Rmd**](Figures.Rmd) is run to create either a bar graph or
boxplot of the global variant-type distributions for each condition in
this Analysis set.
