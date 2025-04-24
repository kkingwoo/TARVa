#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=fp_within
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=00:20:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kkingwoo@charlotte.edu
#SBATCH --error=/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/getRNA_UniqueVarsProgram/fp_within.e
#SBATCH --output=/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/getRNA_UniqueVarsProgram/fp_within.o
#################################################


module load anaconda3
source activate bio_qc

DB_PATH='/projects/kfunk_research2/ChromDBs_AndOther/TARVa.db'
STRING_PATH='/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/DGEFiles_ForStringtie/'
TAB_THREE='/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/getRNA_UniqueVarsProgram/Table3_results.csv'
TPM_FILE='/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/getRNA_UniqueVarsProgram/sig_tpm_diffs.csv'
HGVS_FILE='/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/getRNA_UniqueVarsProgram/Top24_hgvs.txt'
VEP_FILE='/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/getRNA_UniqueVarsProgram/Top24_VEP_out.txt'
DTU_SIG='/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/getRNA_UniqueVarsProgram/tpm_dtu_res_sig.csv'

#python tpm_analysis.py $DB_PATH $STRING_PATH $TAB_THREE $TPM_FILE $HGVS_FILE $VEP_FILE $DTU_SIG

python fpkm_within_condition_analysis.py $DB_PATH $STRING_PATH $TAB_THREE $TPM_FILE $HGVS_FILE $VEP_FILE $DTU_SIG
#If any issues, add cpus-per-task=8 to slurm directives
