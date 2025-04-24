#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=dge1
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kkingwoo@uncc.edu
#SBATCH --error=/scratch/kkingwoo/ADProj/ROSMAP_BulkBrain_2PassMapped/dge1.e
#SBATCH --output=/scratch/kkingwoo/ADProj/ROSMAP_BulkBrain_2PassMapped/dge1.o
#################################################

for dir in `ls -d /projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/Control/*`
do
	theFiles=`ls ${dir}/*fq.gz`
	genome='/users/kkingwoo/FunkLab/STARrefs/Ensembl_hg38Built'
	gtf='/users/kkingwoo/FunkLab/STARrefs/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf'
	pre=${dir##*/}
	prefix=$dir/$pre
	slurm=$prefix'_dge1.slurm'
	job_name=$pre'_dge1'
	err=$prefix'_dge1.e'
	out=$prefix'_dge1.o'
	echo -e "#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name="$job_name "\n#SBATCH --partition=Orion\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=12\n#SBATCH --mem=150gb\n#SBATCH --time=48:00:00\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=kkingwoo@uncc.edu\n#SBATCH --error="$err "\n#SBATCH --output="$out "\n
#################################################\n\nmodule load star\n\nSTAR --runThreadN \$SLURM_NTASKS --genomeDir" $genome "--genomeLoad NoSharedMemory --sjdbGTFfile" $gtf "--readFilesCommand zcat --readFilesIn" $theFiles "--outFileNamePrefix" $prefix "--outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 --chimOutType WithinBAM --sjdbOverhang 139" > $slurm
	sbatch $slurm

done
