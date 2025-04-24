#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=wgs_mub
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --error=wgs_mub.e
#SBATCH --output=wgs_mub.o
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=kkingwoo@uncc.edu
#################################################


for c in `ls -d /nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/*`
do
	for s in `ls -d $c/*`
	do	
		pref=${s##*/}
		pre=${pref%%_Trimmed}
		
		u_job_name=$pre'_unaln'
		u_se=$s/$pre'_unaln.e'
		u_so=$s/$pre'_unaln.o'
		u_slurm=$s/$pre'_unaln.slurm'
		
		echo -e "#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name="$u_job_name "\n#SBATCH --partition=Orion \n#SBATCH --nodes=1 \n#SBATCH --mem=100gb \n#SBATCH --time=24:00:00 \n#SBATCH --error="$u_se "\n#SBATCH --output="$u_so "\n#SBATCH --mail-type=FAIL \n#SBATCH --mail-user=kkingwoo@uncc.edu\n################################################# \n\nmodule load gatk/4.2.0\n\nfor f1 in \`ls -f" $s"/*_R1_tiles.fq.gz\`\ndo\n\tf2=\${f1%%_R1*}'_R2_tiles.fq.gz'\n\tpref=\${f1%%_001_R1*}\n\tpre=\${pref##*/}\n\tunaln=\${pref}'_unaln.bam'\n\thead=\`zcat \$f1 | head -n 1\`\n\thea=\${head% *}\n\tlb="$pre"'-lib'\n\tpu=\${hea#@}\n\n\tgatk FastqToSam -F1 \$f1 -F2 \$f2 -O \$unaln -SM" $pre "-RG \$pre -LB \$lb -PL illumina -PU \$pu\ndone\n\nmodule unload gatk/4.2.0\nmodule load samtools/1.11\n\nout_bam=${s}/${pre}_unaln_mrgd.bam\nfiles=\`ls "$s"/*_unaln.bam\`\nsamtools merge \$out_bam \$files\n\nfor f in \$files\ndo\n\trm \$f\n\ndone" > $u_slurm
		sbatch $u_slurm
	
	done
done
