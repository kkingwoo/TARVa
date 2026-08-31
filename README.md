

**T**ranscriptomic **A**nalysis of **R**NA **Va**riants
  [(TARVa)](#transcriptomic-analysis-of-rna-variants-tarva)                  
  

 ![ ](https://media.giphy.com/media/v1.Y2lkPTc5MGI3NjExOHYxdnZhM3pvdzZxeDY3b3FnbDRxeHJsZHlsNzE1Z3A4MjdqNmZwYyZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/3oKIPnAiaMCws8nOsE/giphy.gif)

# .....UNDER CONSTRUCTION.....                              

## Part 1: Read QC and data processing            


For sorting files by sequencing type and group, and/or to run FastQC and Trimmomatic for read QC. Conda environment dependencies for Part 1 in [bioqc_env.yml](TARVa/bioqc_env.yml). |
### > Sorting input files                  


Run <b>`sort-input`</b> if sequencing files have not yet been sorted and organized by group/condition, WGS/RNAseq.               
            
--<b>Usage:</b> `tarva sort-input --metadata <metadata.csv> --input <input_dir> --output <sorted_dir>`


### > Read QC using Trimmomatic and FastQC   


Run <b>`qc`</b> to apply FastQC and Trimmomatic to the data.          


--<b>Usage:</b>```tarva qc --metadata <metadata.csv> --input <fastq_dir> --output <output_dir>```             

For full set of options ```tarva qc --help```                 

When running this step, QC output file locations will be written to  <b>*tarva_qc_output_location.txt* </b> in the path. 

<details>                
<summary><b>Click to see details about metadata file for input</b></summary>  


The modules used in this part require a *metadata.csv* file with the following format:            


| group | individualID | WGS | RNA | wgs_barcode | rna_barcode | 
|---|---|---|---|---|---|---|
| male | m1 | ERR2814681 | ERR2812923 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa |
| male | m2 | ERR2814680 | ERR2812921 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa |
| male | m3 | ERR2814679 | ERR2812918 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa |
| female | f1 | ERR2814677 | ERR2812930 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa |
| female | f2 | ERR2814678 | ERR2812929 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa |
| female | f3 | ERR2814676 | ERR2812927 | TruSeq3-PE-2.fa | TruSeq3-PE-2.fa |                   

1.) Exactly two distinct groups/conditions should be present in the 'group' column.                
2.) individualID is the UUID used for each individual.                  
3.) Only the unique prefixes from raw fastq file identifiers (columns 'WGS' and 'RNA') should be the values, ignoring any specific string patterns that identify paired-end reads, (e.g., "_1","_2", "R1", "R2").                               
4.) Barcode/adapter files for use with Trimmomatic are required for QC steps.                 

</details> 

### > Read alignment                

WGS alignment is carried out using BWA-mem. RNAseq alignment is run twice:  once for Stringtie, and again in two-pass mode, for variant calling. 

### > Alignment QC                


### > Variant calling and RNA-> WGS check                 


## Part 2: Local DB Build                

For sorting files by sequencing type and group, and/or to run FastQC and Trimmomatic for read QC. Conda environment dependencies for Part 2 in [tarva_env.yml](TARVa/tarva_env.yml).            


## Part 3: Data analysis                  
<details>                
<summary><b>Click to see details about sample info file for input</b></summary> 

A *sample_info.csv* file is required for this part. 'individualID' field is required. Any variables to be tested in relation to the RNA variants identified should each have their own column. See sample file below, where 'Age' is the only variable that will be assessed.         


| individualID | Age(mos) |                  
|---|---|                 
| m1 | 2 |
| m2 | 10 |
| m3 | 3 |
| f1 | 1 |
| f2 | 8 |
| f3 | 5 |

</details> 
--------------------------------------------------------------------
--------------------------------------------------------------------
> [ * * * <b>CLI WIP</b> * * * ]               

<b>Pipelines to test and incorporate:</b>                 

- [X] Data QC                  
- [] Data processing                  
- [] DB build                 
- [] Data analysis                          
- [] Finishing touches

---------------------------------------------------------------
---------------------------------------------------------------

### *We will return after the job is done. Stay tuned*                  


<u>AI assistance with code</u>                    
                
**Original code and logic:**
[![Quillx](https://raw.githubusercontent.com/qainsights/Quillx/main/badges/quillx-1.svg)](https://github.com/qainsights/Quillx)

**Original code and logic --> CLI:**
[![Quillx](https://raw.githubusercontent.com/qainsights/Quillx/main/badges/quillx-3.svg)](https://github.com/qainsights/Quillx)                 

### References                


                    
                    

[License](LICENSE)

