#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=final_outs
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --error=final_outs.e
#SBATCH --output=final_outs.o
#################################################

module load anaconda3
source activate tarva

DB_PATH='TARVa.db'
GENES_FILE='adj_all_proportions_pvals_0_075_filter.csv'
A1_FILE='alt1_per_position_gene.csv'
A2_FILE='alt2_per_position.csv'
REDI_DB='TABLE1_hg38.txt'


python anet.py $DB_PATH $GENES_FILE $A1_FILE $A2_FILE $REDI_DB 
