#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=mm_star
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --error="star_par.e"
#SBATCH --output="star_par.o"
#################################################

module load star

GTF="gencode.v42.annotation.gff3"
GENOME="final_run_gencode_hg38Built"
DICT="dict_dumpfi.json"
FASTA="GRCh38.primary_assembly.genome.fa"
GATK_SLURM="full_gatk.slurm"
DUMPOUT="gatk_tupes.json"

python star.py "$GENOME" "$GTF" "$DICT" "$SLURM_CPUS_PER_TASK" "$GATK_SLURM" "$DUMPOUT"

