#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=dge2_3
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#SBATCH --error=dge2_3e
#SBATCH --output=dge2_3.o
#################################################

for in_bam in `ls -d /path/to/RNAseq/*/*/*_differential_Aligned.sortedByCoord.out.bam`
do	
	gtf='All_mrgd.gtf'
	#gtf='/users/kkingwoo/FunkLab/STARrefs/'
	pr=${in_bam%/*}
	#rm $pr/*tab
	#rm $pr/*gtf
	pre=${pr##*/}
	pref=${pr%%differential_Aligned.sortedByCoord.out.bam}/$pre
	slurm=$pref'_dge2_3_fasta.slurm'
	job_name=$pre'_dge2_3_fasta'
	err=$pref'_dge2_3_fasta.e'
	out=$pref'_dge2_3_fasta.o'
	out_trans=$pref'_round3_fasta.gtf'
	tab=$pref'_round_fasta.tab'
	full_cov=$pref'_r3_cov_fasta.gtf'
	echo -e "#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name="$job_name "\n#SBATCH --partition=\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=12\n#SBATCH --mem=64gb\n#SBATCH --time=12:00:00\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=\n#SBATCH --error="$err "\n#SBATCH --output="$out "\n#################################################\n\nstringtie $in_bam -o $out_trans --rf -m 50 -t -p \$SLURM_NTASKS -G $gtf -A $tab -C $full_cov -B -e" > $slurm
	sbatch $slurm

done
