#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=dge2_3
#SBATCH --partition=Nebula
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=150gb
#SBATCH --time=5:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kkingwoo@uncc.edu
#SBATCH --error=/scratch/kkingwoo/ADProj/ROSMAP_BulkBrain_2PassMapped/dge2_3e
#SBATCH --output=/scratch/kkingwoo/ADProj/ROSMAP_BulkBrain_2PassMapped/dge2_3.o
#################################################

python3 DGE_Pipe2_Round2.py 
