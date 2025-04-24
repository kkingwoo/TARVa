import csv
import os


overlapped = open('/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/RNASeqFiltered_VCFs_withTBI/OverlappedSamples.txt', 'r')
overlap = csv.reader(overlapped, delimiter='\t')
overlapList = []

for o in overlap:
    subdir=o[0]
    overlapList.append(subdir)

intervals = '/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/somatic-hg38_CNV_and_centromere_blacklist.hg38liftover.list'
ref = '/users/kkingwoo/FunkLab/BWA_Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
dic = '/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict'
gr = '/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/somatic-hg38_af-only-gnomad.hg38.vcf'
pons='/nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/pons.vcf.gz'

main_dir = '/nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/AD/'

for a in os.listdir(main_dir):
    ai = a.split('_')[0]
    if ai in overlapList:
        a_path = main_dir+a
        slurm=a_path+'/'+ai+'_Pons_toAD.slurm'
        job_name=ai+'_Pons_toAD'
        out=a_path+'/'+ai+'_Pons_toAD.o'
        err=a_path+'/'+ai+'_Pons_toAD.e'
        in_bam=a_path+'/'+ai+'_recal.bam'
        out_vcf=a_path+'/'+ai+'_pons.vcf.gz'


        os.system("echo -e '#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name={0}\n#SBATCH --partition=Orion\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=4\n#SBATCH --mem=150gb\n#SBATCH --time=1-00:00:00\n#SBATCH --error={1}\n#SBATCH --output={2}\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=kkingwoo@uncc.edu\n#################################################\n\nmodule load gatk/4.2.0\n\ngatk Mutect2 -R {3} -I {4} -tumor {5} --germline-resource {6} --panel-of-normals {7} --max-mnp-distance 0 -O {8}' > {9}".format(job_name,err,out,ref,in_bam,ai,gr,pons,out_vcf,slurm))
        os.system("sbatch {0}".format(slurm))
