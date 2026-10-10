

### **T**ranscriptomic **A**nalysis of **R**NA **Va**riants [(TARVa)](#transcriptomic-analysis-of-rna-variants-tarva)  
  
## Project Overview                                   
  [**TARVa**](#transcript-analysis-of-rna-variants-tarva) is a tool which enables analysis of variant calls from raw RNA sequencing data, taking  into consideration both depth of read and length of transcript, and scans for all variant call-types in a strand-aware fashion. 
  Starting with a global overview, the pipeline moves through non-hierarchical stratification of the data, ending with local, per-position comparisons across tissue types and conditions. 
  The  workflow is set up for use with paired-sample raw RNAseq and whole genome sequencing files (fastq) from more than one condition and/or tissue from diploid organisms,
  and covers data preprocessing in addition to the main TARVa analysis pipeline.                
  __NOTE:__ Parts or all of the workflow, algorithm(s), and other components of this project are patent pending.  
  
  The workflow contains three main pipelines to be run consecutively:               
  1. [Read QC and Data Processing](#part-one)                 
  2. [Local Database Build](#part-two)                
  3. [Data Analysis](#part-three)                   
  
  The project was originally created specifically for analysis of bulk brain tissue samples 
  from Alzheimer's and cognitively normal aged individuals. The fourth section, [Optional Modules](#part-optional) 
  includes analyses included in the manuscript but that may not be applicable to other tissue types
  or organisms and thus not critical for the TARVa workflow.                      
  
  Due to resource requirements, it is recommended to use the tool's `--hpc` mode whenever possible.

## Tools, Environments and Dependencies
Python 3 is a requirement for this workflow, and three non-interchangeable conda environments are used: [__*bio_qc*__](bioqc_env.yml), [__*stringtie*__](stringtie_env.yml), [__*tarva*__](tarva_env.yml)                                                
Installation of __star/2.7.11b__, __bwa/0.7.17__, and __R/4.5.1__ are also needed.            

__Python modules__              
  - pandas            
  - numpy           
  - scipy           
  - statsmodels                 
  - matplotlib            
  

 ![ ](https://media.giphy.com/media/v1.Y2lkPTc5MGI3NjExOHYxdnZhM3pvdzZxeDY3b3FnbDRxeHJsZHlsNzE1Z3A4MjdqNmZwYyZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/3oKIPnAiaMCws8nOsE/giphy.gif)

# .....UNDER CONSTRUCTION.....                             

## <b>1. <a id="part-one"></a>Read QC and Data Processing Pipeline </b>                  

This pipeline accepts paired-end sequencing FASTQ files and returns RNA-specific variant calls. 
Reads are sorted by sequencing type and group, then quality-trimmed and aligned.  __BWA-mem__ [[1]](#ref-bwa) is used for WGS data and __STAR__ [[2]](#ref-star) for RNA-seq.                
The RNA-seq alignments are also assembled with __StringTie__ [[3]](#ref-stringtie) which allows evaluation of gene expression levels and enables creation of reference files containing both known and novel transcripts.             

Data for both sequencing types are run through __GATK__ [[4]](#ref-gatk) short-variant discovery pipelines. Two filtered VCF files per individual are created, one for each
sequencing type.   

Activate the [bio_qc](bioqc_env.yml) environment for the first part of this pipeline. The [tarva](tarva_env.yml) environment should be activated for the last step.               


Input Files:                    

  - Reference genome FASTA                        
  - Reference genome GTF and GFF3                 
  - Short-read RNA-seq and WGS-seq FASTQ sequencing files for each individual from two different conditions (minimum # samples recommended n=30 for each condition)               
  - Metadata CSV file.  See [NOTE](#note-metadata) below for more information.                    
  

### <b><a id="note-metadata"></a>NOTE:</b> metadata input                     

<details>                
<summary><b>Click to see details about metadata file for input</b></summary>  


The modules used in this part require a *metadata.csv* file. 'individualID' field is required. 
Any variables to be tested in relation to the RNA variants, in addition to the 'group' field, should each have their own column. 
See sample file below, where 'Age' is the meta-variable that will be assessed. __NOTE__: The data analysis pipeline has sufficient power when sample sizes are n=29 or larger.                        


| group | individualID | WGS | RNA | wgs_barcode | rna_barcode | age(mos) |
|---|---|---|---|---|---|---|
| control | c1 | R281468 | R281293 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa | 2 |   
| control | c2 | R281460 | R281291 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa | 10 |   
| control | c3 | R281467 | R281298 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa | 3 |   
| treatment | t1 | R281462 | R281290 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa | 1 |   
| treatment | t2 | R281465 | R281299 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa | 8 |   
| treatment | t3 | R281469 | R281297 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa | 5 |                   

1.) Exactly two distinct groups/conditions/tissues should be present in the 'group' column.                
2.) individualID is the UUID used for each individual.                  
3.) Only the unique prefixes from raw fastq file identifiers (columns 'WGS' and 'RNA') should be the values, ignoring any specific string patterns that identify paired-end reads, (e.g., "_1","_2", "R1", "R2").                               
4.) Barcode/adapter files for use with Trimmomatic are required for QC steps.                 

</details>                

</blank><blank>                     

### <b>NOTE:</b> project directory tree                    

<details>             

<summary><strong>Expected Directory Structure (click to expand)</strong></summary>

`<group>` and `<id>` are values from the metadata `group` and `individualID` columns; `<adapter>` is the adapter FASTA's filename without `.fa`.

```
<project_root>/
├── metadata.csv                            INPUT
├── RawData/                                INPUT -- any layout, searched recursively
│   └── **/<prefix>*.fastq.gz
├── Refs/                                   INPUT
│   ├── <genome>.fna                          + .fai and .dict from prep-reference
│   ├── <genome>.gtf                          STAR annotation
│   ├── <genome>.gff                          StringTie / gffcompare annotation
│   ├── star_<genome>/                        prep-reference
│   ├── star_<genome>_combined/               prep-reference, on combined_reference.fa
│   └── KnownSites/
│       └── <known_sites>.vcf.gz
│
├── Sorted/                                 sort-input
│   ├── WGS/<group>/<adapter>/*.fastq.gz
│   └── RNA/<group>/<adapter>/*.fastq.gz
│
├── QC/                                     qc
│   ├── fastqc/WGS|RNA/<group>/<adapter>/
│   ├── trimmed/WGS|RNA/<group>/<adapter>/*.paired.fastq.gz
│   └── *_summary.csv, read_lengths.csv
│
├── STAR_FirstPass/                         stringtiepipe align
│   ├── <group>/<id>/<id>_firstpass.bam
│   └── *_summary.csv
│
├── FirstPassAssembly/                      stringtiepipe assemble
│   ├── merged.gtf
│   ├── <group>/<id>/<id>_round{1,3}.gtf
│   ├── <group>/<id>/*.ctab
│   └── *_summary.csv
│
├── NovelTranscripts/                       stringtiepipe novel-transcripts
│   ├── gffcompare/
│   ├── novel.gtf
│   ├── novel_transcripts.fa
│   └── combined_reference.fa                 + .fai and .dict
│
├── WGS_Aligned/                            wgseqpipe align
│   ├── <group>/<id>/<id>_wgs.bam
│   └── *_summary.csv
│
├── RNA_TwoPass/                            rnaseqpipe align
│   ├── <group>/<id>/<id>_twopass.bam
│   └── *_summary.csv
│
├── WGS_Preprocessed/                       wgseqpipe preprocess
│   ├── <group>/<id>/<id>_bqsr.bam
│   └── *_summary.csv
│
├── RNA_Preprocessed/                       rnaseqpipe preprocess
│   ├── <group>/<id>/<id>_bqsr.bam
│   └── *_summary.csv
│
├── WGS_PON/                                wgseqpipe pon
│   ├── calls/<id>_pon.vcf.gz
│   ├── pon_db/
│   ├── pon.vcf.gz
│   └── *_summary.csv
│
├── WGS_Variants/                           wgseqpipe variants
│   ├── <group>/<id>/<id>_contaminates_filtered.vcf.gz
│   └── *_summary.csv
│
├── RNA_Variants/                           rnaseqpipe variants
│   ├── <group>/<id>/<id>_varcall_flt.vcf
│   └── *_summary.csv
│
└── RNA_Unique/                             analysispipe rna-unique
    ├── <group>/<id>/<id>_rna_unique.vcf
    └── *_summary.csv
```

Every per-sample step also writes `logs/` and, with `--array`, `summary_parts/`.

</details>                

</blank><blank>     

### a. <u>Sort input files and prepare reference files</u>                 


- Sort and organize input sequencing files if they have not yet been sorted and organized by group/condition, WGS/RNAseq.                  

<b>Usage:</b> `tarva sort-input --help`      


- Index and prepare reference files for WGS and RNA-seq alignment tool compatibility.  Accepts compressed and uncompressed formats.           


<b>Usage:</b> `tarva prep-reference --help`                 



### b. <u>Read QC</u>   


- Apply FastQC and Trimmomatic to the data, remove low-quality tiles as needed, and return trimmed fastq files.                 
  This step can be repeated as many times as needed. Read lengths are assessed and metric summaries recorded for each job.                


<b>Usage:</b>```tarva qc --help```                             



### c. <u>Stringtie</u>            

  - <b>Step 1:</b> First pass alignment with STAR.                  
  
  <b>Usage:</b> ```tarva stringtiepipe align --help ```                   
  
  - <b>Step 2:</b> Assemble transcripts with Stringtie using STAR-aligned reads from Step 1.                  
    `--strandedness` is required, see [[5]](#stringtie-manual) for Stringtie format if data are paired-end.                   
  
  <b>Usage:</b> ```tarva stringtiepipe assemble --help ```                  
  
  - <b><a id="step-three"></a>Step 3:</b> Identify novel transcripts and merge Stringtie-created reference files with original reference files.                
    Output files are used for subsequent data processing steps.                 
  
  <b>Usage:</b> ```tarva stringtiepipe novel-transcripts --help ```                 
  
  
### d. <u>WGS</u>                 

  - <b>Step 4:</b> Align WGS reads using BWA-mem.                     
  
  <b>Usage:</b> ```tarva wgseqpipe align --help ```                      
  
  - <b>Step 5:</b> QC alignments with GATK pipeline.            
  
  <b>Usage:</b> ```tarva wgseqpipe preprocess --help ```                    
  
  - <b>Step 6:</b> Create PONs from control group (Optional).
    Provides the option to create a PONs with the study's control group for use in [Step 7](#step-seven)                  
  
  <b>Usage:</b> ```tarva wgseqpipe pon --help ```               
  
  - <b><a id="step-seven"></a>Step 7:</b> WGS variant calling. Output VCF files.               
  Uses Mutect2 with optional --germline-resource and/or --pon panel of normals (PONs).                  
  IMPORTANT NOTE: If one or either option is used in this step, masking must be done for the matching sites where 
  RNA variants were identified in [Step 11](#step-eleven). Application of PONs to WGS-only causes ambiguity 
  around origin of variant and may inflate RNA-specific assignment of variant calls.                
  
  <b>Usage:</b> ```tarva wgseqpipe variants --help ```                      
  

### e. <u>RNA-seq</u>                   

  - <b>Step 8:</b> Build index with `prep-reference` for combined reference created by Stringtie in [Step 3](#step-three).                 

  - <b>Step 9:</b> Align RNA-seq reads using BWA-mem.                     
  
  <b>Usage:</b> ```tarva rnaseqpipe align --help ```                      
  
  - <b>Step 10:</b> QC alignments with GATK short-variant pipeline.                      
  
  <b>Usage:</b> ```tarva rnaseqpipe preprocess --help ```                    
  
  - <b><a id="step-eleven"></a>Step 11:</b> WGS variant calling. Output VCF files.                     
  
  <b>Usage:</b> ```tarva rnaseqpipe variants --help ```

### f. <u>Identify RNA-specific variants</u>                 

  - <b>Step 11:</b> Match RNA and WGS variant call sites for each individual.               
    If PON or germline resource were used during WGS `--align`, those sites should also be masked in RNA
    to reduce false positives.                
    
  <b>Usage:</b> ```tarva analysispipe rna-unique --help ```
  

            


## <b>2. <a id="part-two"></a>Local DB Build</b>                



## <b>3. <a id="part-three"></a>Data Analysis</b>   

## <b>4. <a id="part-optional"></a>Optional Modules</b>             




---------------------------------------------------------------
---------------------------------------------------------------                

### [ * * * <b>CLI WIP</b> * * * ]               

<b>Pipelines tested and incorporated:</b>                 

- [X] Data QC (10/9/2026)                 
- [X] Data processing (10/9/2026)                 
- [] DB build 
- [] Data analysis 
- [] Other              


---------------------------------------------------------------
---------------------------------------------------------------

### *We will return after the job is done. Stay tuned*                  


<u>AI assistance with code</u>                    
                
**Original code and logic:**
[![Quillx](https://raw.githubusercontent.com/qainsights/Quillx/main/badges/quillx-1.svg)](https://github.com/qainsights/Quillx)

**Original code and logic converted to CLI:**
[![Quillx](https://raw.githubusercontent.com/qainsights/Quillx/main/badges/quillx-3.svg)](https://github.com/qainsights/Quillx)                 

### References          
1. <a id="ref-bwa"></a>Li, Heng, and Richard Durbin. “Fast and Accurate Short Read Alignment with Burrows–Wheeler Transform.” Bioinformatics 25, no. 14 (2009): 1754–60. https://doi.org/10.1093/bioinformatics/btp324.            

2. <a id="ref-star"></a>Dobin, Alexander, Carrie A. Davis, Felix Schlesinger, et al. “STAR: Ultrafast Universal RNA-Seq Aligner.” Bioinformatics 29, no. 1 (2013): 15–21. https://doi.org/10.1093/bioinformatics/bts635.            

3. <a id="ref-stringtie"></a>Shumate, Alaina, Brandon Wong, Geo Pertea, and Mihaela Pertea. “Improved Transcriptome Assembly Using a Hybrid of Long and Short Reads with StringTie.” PLOS Computational Biology 18, no. 6 (2022): e1009730. https://doi.org/10.1371/journal.pcbi.1009730.            

4. <a id="ref-gatk"></a>Brouard, Jean-Simon, Flavio Schenkel, Andrew Marete, and Nathalie Bissonnette. “The GATK Joint Genotyping Workflow Is Appropriate for Calling Variants in RNA-Seq Experiments.” Journal of Animal Science and Biotechnology 10, no. 1 (2019): 44. https://doi.org/10.1186/s40104-019-0359-0.            

5. <a id="stringtie-manual"></a><https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual>



                    
                    

[License](LICENSE)

