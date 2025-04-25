#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=pileups
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --error=pileups.e
#SBATCH --output=pileups.o
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=
#################################################


python3 Step9_GetPileupSummaries_Control.py

python3 Step9_GetPileupSummaries_AD.py

