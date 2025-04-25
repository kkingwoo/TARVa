#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=fp_within
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=00:20:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --error=fp_within.e
#SBATCH --output=fp_within.o
#################################################


module load anaconda3
source activate bio_qc

DB_PATH='TARVa.db'
STRING_PATH='DGEFiles_ForStringtie/'
TAB_THREE='Table3_results.csv'
TPM_FILE='sig_tpm_diffs.csv'
HGVS_FILE='Top24_hgvs.txt'
VEP_FILE='Top24_VEP_out.txt'
DTU_SIG='tpm_dtu_res_sig.csv'

#python tpm_analysis.py $DB_PATH $STRING_PATH $TAB_THREE $TPM_FILE $HGVS_FILE $VEP_FILE $DTU_SIG

python fpkm_within_condition_analysis.py $DB_PATH $STRING_PATH $TAB_THREE $TPM_FILE $HGVS_FILE $VEP_FILE $DTU_SIG
#If any issues, add cpus-per-task=8 to slurm directives
