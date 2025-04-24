#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=create_pon
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --error=create_pon.e
#SBATCH --output=create_pon.o
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=kkingwoo@uncc.edu
#################################################

#module load gatk/4.2.0

#gatk IndexFeatureFile -I /users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/somatic-hg38_af-only-gnomad.hg38.vcf 

#python3 Step7.py

#ref='/users/kkingwoo/FunkLab/BWA_Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

python3 Step7_CreatePons.py

