#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=full_gatk
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=128gb
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --error="/path/to/full_gatk.e"
#SBATCH --output="/path/to/full_gatk.o"
#################################################

module load anaconda3
source activate bio_qc
module load R/4.3.3
export R_LIBS_USER=Rlibs

DUMPED="gat_tupes.json"
FA="GRCh38.primary_assembly.genome.fa"
DS_VAR="monothreecon_bulkbrainmci"
KNOWN="resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"

python full_gatk.py "$DUMPED" "$FA" "$DS_VAR" "$KNOWN"

