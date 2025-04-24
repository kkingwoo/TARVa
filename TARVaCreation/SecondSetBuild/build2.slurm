#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=combo
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --mem=8gb
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=
#SBATCH --error=combo.e
#SBATCH --output=combo.o
#################################################

## for running all Build_SecondSet.py functions, add --cpus-per-task=16 and --mem=128gb ##

module load anaconda3
source activate tarva

DB_PATH='TARVa.db'
MONO_PATH='/path/to/monocyte/vcfs/'
BB_PATH='/path/to/brain/vcfs/'
CLIN='ROSMAP_clinical.csv'
OUT='secondset_out.csv'
BB_FIRST='Table3_results.csv'
MOD_DUMP='top78mods.json'
PROF_MODS_IN='for_profile_input_mods.csv'
PROF_MODS_OUT='top24_AD_profile.csv'
JSON1='make_profs.json'
JSON2='samp_dict.json'
JSON3='prof_deets.json'
TOPVARS='78Vars.csv'

python Build_SecondSet.py $DB_PATH $MONO_PATH $BB_PATH $CLIN $OUT $BB_FIRST $MOD_DUMP $PROF_MODS_IN $PROF_MODS_OUT $JSON1 $JSON2 $JSON3 $TOPVARS
