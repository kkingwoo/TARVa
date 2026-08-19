#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=final_outs
#SBATCH --partition=Nebula
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --error=final_outs.e
#SBATCH --output=final_outs.o
#################################################

module load anaconda3
source activate tarva

DB_PATH='TARVa.db'
METADAT='metada.csv'


python anet.py $DB_PATH $METADAT 
