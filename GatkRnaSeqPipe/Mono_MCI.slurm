#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=mono_mci
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#SBATCH --error=mono_mci.e
#SBATCH --output=mono_mci.o
#################################################
module load anaconda3
source activate bio_qc

DB_PATH='TARVa.db'
MONO_MAN='SYNAPSE_METADATA_MANIFEST_Mono.tsv' #list of sample data retrieved from ROSMAP download of blood/monocyte data
CLIN_DAT='ROSMAP_clinical.csv' #Clinical information about samples (in this case, specific to ROSMAP)
BB_MAN='SYNAPSE_METADATA_MANIFEST_BB.tsv' #list of sample data retrieved from ROSMAP download of bulk brain data
BB_FAST='/path/to/tissue1_fastqs/' 
MONO_FAST='/path/to/tissue2_fastqs/'
GTF='gencode.v42.annotation.gff3'
GENOME='final_run_gencode_hg38Built'
SLURM='STAR_array.slurm'
DUMPS='dict_dumpfi.json'
GATK_SLURM='full_gatk.slurm'
DUMPOUT='gat_tupes.json'


python match_monocytes.py $DB_PATH $MONO_MAN $CLIN_DAT $BB_MAN $BB_FAST $MONO_FAST $GTF $GENOME $SLURM $DUMPS $GATK_SLURM $DUMPOUT
