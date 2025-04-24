#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=GP1
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3-12:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=
#SBATCH --error=gp1.e
#SBATCH --output=gp1.o
#################################################

module load anaconda3
source activate tarva

GTF_FILE_PATH='All_mrgd.gtf'
MERGED_FASTA='gencode_stringtie_merged_output.fasta'
DB_PATH='TARVa.db'
GENES_FILE='GroupedByGenesRawSigPvals.csv'
SINGLE_ISO_FILE='single_iso_genes.csv'
WGS_VCF_PATH='WGS_VCFs/'
RNA_VCF_PATH='RNA_VCFs/'
STRING_PATH='DGEFiles_ForStringtie/'
COUNTS_FILE='per_position_gene.csv'

python BuildTarvaDBs.py $GTF_FILE_PATH $MERGED_FASTA $DB_PATH $GENES_FILE $SINGLE_ISO_FILE $WGS_VCF_PATH $RNA_VCF_PATH $STRING_PATH $COUNTS_FILE 
