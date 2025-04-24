#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=md_ss
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --error=md_ss.e
#SBATCH --output=md_ss.o
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=kkingwoo@uncc.edu
#################################################

python3 Step5_MD_SS.py
