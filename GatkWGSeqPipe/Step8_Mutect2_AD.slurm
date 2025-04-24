#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=apply_pon
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --error=apply_pon.e
#SBATCH --output=apply_pon.o
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=kkingwoo@uncc.edu
#################################################


python3 Step8_ApplyPonsToAD.py

