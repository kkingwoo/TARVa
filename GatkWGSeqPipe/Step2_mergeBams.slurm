#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=wgs_smb
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --error=wgs_smb.e
#SBATCH --output=wgs_smb.o
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=
#################################################

module load samtools/1.11

for c in `ls -d /path/to/ROSMAP_WGS/TrimmedFASTQFiles/*`
do
	for s in `ls -d $c/*`
	do
		for bd in `ls -d $s/*/`
		do
			samp=${bd##*_Trimmed/}
			sa=${samp%%_*}
			files=`ls $bd*_srt.bam`
			pre=$sa'_merged.bam'
			merged_bam=$bd$pre
			job_name=$sa'_merge'
			se=$bd$sa'_merge.e'
			so=$bd$sa'_merge.o'
			slurm=$bd$sa'_merge.slurm'

				
			echo -e "#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name="$job_name "\n#SBATCH --partition= \n#SBATCH --nodes=1 \n#SBATCH --mem=64gb \n#SBATCH --time=24:00:00 \n#SBATCH --error="$se "\n#SBATCH --output="$so "\n#SBATCH --mail-type=FAIL \n#SBATCH --mail-user=\n################################################# \n\nmodule load samtools/1.11\n\nsamtools merge "$merged_bam $files "\n\nfor f in "$files "\n\tdo\n\t\trm \$f\n\ndone" > $slurm
			sbatch $slurm
		
	
		done
	done
done
