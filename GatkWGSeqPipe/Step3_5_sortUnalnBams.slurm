#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=wgs_mub
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --error=wgs_mub.e
#SBATCH --output=wgs_mub.o
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=
#################################################


for c in `ls -d /path/to/ROSMAP_WGS/TrimmedFASTQFiles/*/*`
do
	
	pref=${c##*/}
	pre=${pref%%_Trimmed}
	ss_job_name=$pre'_ss'
	ss_se=$c/$pre'_ss.o'
	ss_so=$c/$pre'_ss.e'
	ss_slurm=$c/$pre'_ss.slurm'
	
	echo -e "#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name="$ss_job_name "\n#SBATCH --partition= \n#SBATCH --nodes=1 \n#SBATCH --ntasks-per-node=16 \n#SBATCH --mem=16gb \n#SBATCH --time=24:00:00 \n#SBATCH --error="$ss_se "\n#SBATCH --output="$ss_so "\n#SBATCH --mail-type=FAIL \n#SBATCH --mail-user=\n#################################################\n\nmodule load gatk/4.2.0\n\nfor umb in \`ls -f" $c"/*_unaln_mrgd.bam\`\ndo\n\tfile_name=\${umb##*/}\n\tdir_name=\${umb%%\$file_name}\n\ttmp=\${dir_name}'TmpDir'\n\tmkdir \$tmp\n\texport _JAVA_OPTIONS=-Djava.io.tmpdir=\$tmp\n\to_umb=\${umb%%.bam}'_srtd.bam'\n\tgatk SortSam -I \$umb -O \$o_umb -SO queryname --TMP_DIR \$tmp\nrm \$umb\ndone" > $ss_slurm
	sbatch $ss_slurm
	
	
done
