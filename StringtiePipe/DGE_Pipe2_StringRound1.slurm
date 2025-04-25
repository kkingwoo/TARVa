#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=dge1
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=
#SBATCH --error=dge1.e
#SBATCH --output=dge1.o
#################################################

for in_bam in `ls -d /path/to/RNAseq/*/*/*_differential_Aligned.sortedByCoord.out.bam`
do
	gtf='gencode.v42.annotation.gff3'
	pr=${in_bam%/*}
	pre=${pr##*/}
	pref=${pr%%Aligned.sortedByCoord.out.bam}/$pre
	slurm=$pref'_dge2.slurm'
	job_name=$pre'_dge2'
	err=$pref'_dge2.e'
	out=$pref'_dge2.o'
	#out_bam=$pref'_star_srtd.bam'
	out_trans=$pref'_round1.gtf'
	tab=$pref'_round1.tab'
	full_cov=$pref'_cov_round1.gtf'

	echo -e "#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name="$job_name "\n#SBATCH --partition=\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=12\n#SBATCH --mem=64gb\n#SBATCH --time=24:00:00\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=\n#SBATCH --error="$err "\n#SBATCH --output="$out "\n#################################################\n\nrm *ctab\n\nstringtie $in_bam -o $out_trans --rf -t -p \$SLURM_NTASKS -G $gtf -A $tab -C $full_cov -B " > $slurm
	sbatch $slurm

done
