import os
import csv

known='/users/kkingwoo/FunkLab/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf'
ref = '/users/kkingwoo/FunkLab/BWA_Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

overlapped = open('/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/RNASeqFiltered_VCFs_withTBI/OverlappedSamples.txt', 'r')
overlap = csv.reader(overlapped, delimiter='\t')
overlapList = []

pons = open('/nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/pons.list', 'w')
pons_write = csv.writer(pons, delimiter='\t')

gr = '/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/somatic-hg38_af-only-gnomad.hg38.vcf'

for o in overlap:
    subdir=o[0]
    overlapList.append(subdir)

main_dir = '/nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/Control/'

for c in os.listdir(main_dir):
    ci = c.split('_')[0]
    if ci in overlapList:
        c_path = main_dir+c
        job_name=ci+'_create_pon'
        error=c_path+'/'+ci+'_create_pon.e'
        out=c_path+'/'+ci+'_create_pon.o'
        slurm_script=c_path+'/'+ci+'_create_pon.slurm'
        in_bam=c_path+'/'+ci+'_recal.bam'
        out_vcf=c_path+'/'+ci+'_pon.vcf.gz'
        pons_write.writerow([out_vcf])
        print(c_path)
        os.system("echo -e '#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name={0}\n#SBATCH --partition=Orion\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=4\n#SBATCH --mem=150gb\n#SBATCH --time=1-00:00:00\n#SBATCH --error={1}\n#SBATCH --output={2}\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=kkingwoo@uncc.edu\n#################################################\n\nmodule load gatk/4.2.0\n\ngatk Mutect2 -R {3} -I {4} -tumor {5} --germline-resource {6}' --max-mnp-distance 0 -O {7} > {8}".format(job_name,error,out,ref,in_bam,ci,gr,out_vcf,slurm_script))

        os.system("sbatch {0}".format(slurm_script))

overlapped.close()
