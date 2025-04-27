#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=dsan
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --error=dsan.e
#SBATCH --output=dsan.o
#################################################

module load anaconda3
source activate tarva

DB_PATH='TARVa.db'
AD_ONLY='AD_only_For_Enrichr.txt'
AD_ONLY_PATTS='primarily_ad_patterns.csv'
KNOWN_V_NOVEL='primarily_ad_known_and_novels.csv'
FILTERED_LOCAL='local_level_highly_sig_type_diffs_per_gene.csv'

python downstream_analysis.py $DB_PATH  $AD_ONLY $AD_ONLY_PATTS $KNOWN_V_NOVEL $FILTERED_LOCAL
