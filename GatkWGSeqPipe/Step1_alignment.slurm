#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=wgs_bma_t
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --error=wgs_bma_t.e
#SBATCH --output=wgs_bma_t.o
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=kkingwoo@uncc.edu
#################################################

ref=/users/kkingwoo/FunkLab/BWA_Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz



for c in `ls -d /nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/*`
do
	for s in `ls -d $c/*`
	do
		
		sa=${s##*/}
		sample=${sa%%_*}
		sort_dir=$s/$sample'_SortedBams'
		mkdir -p $sort_dir
		slurm=$s/$sample'_dirSub.slurm'
		job_name=$sample'_dirSub'
		so=$s/$sample'_dirSub.o'
		se=$s/$sample'_dirSub.e'
	 		
		echo -e "#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name="$job_name "\n#SBATCH --partition=Orion \n#SBATCH --nodes=1 \n#SBATCH --ntasks-per-node=16 \n#SBATCH --mem=16gb \n#SBATCH --time=24:00:00 \n#SBATCH --error="$se "\n#SBATCH --output="$so "\n#SBATCH --mail-type=FAIL \n#SBATCH --mail-user=kkingwoo@uncc.edu\n################################################# \n\nmodule load bwa\n\nmodule load samtools\n\nfor f1 in \`ls -f" $s"/R*R1*fq.gz\`\ndo\n\tf2=\${f1%%R1_tiles*}R2_tiles.fq.gz\n\tprefix=\${f1##*/}\n\tpre=\${prefix%%_001_*}\n\tout_sam="$s"/\$pre.sam\n\tout_bam="$s"/\$pre.bam\n\tsrt_bam="$sort_dir"/\${pre}_srt.bam\n\n\tbwa mem -t \$SLURM_NTASKS $ref \$f1 \$f2 > \$out_sam\n\tsamtools view -S -b \$out_sam > \$out_bam\n\trm \$out_sam\n\tsamtools sort \$out_bam -o \$srt_bam\n\trm \$out_bam\ndone" > $slurm
		sbatch $slurm
		
	done
done
