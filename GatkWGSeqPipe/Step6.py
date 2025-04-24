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
import os
import csv

known='/users/kkingwoo/FunkLab/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf'
ref = '/users/kkingwoo/FunkLab/BWA_Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'

overlapped = open('/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/RNASeqFiltered_VCFs_withTBI/OverlappedSamples.txt', 'r')
overlap = csv.reader(overlapped, delimiter='\t')
overlapList = []

for o in overlap:
    subdir=o[0]
    overlapList.append(subdir)

cList = ["AD", "Control"]

main_dir = '/nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/'

for c in os.listdir(main_dir):
    if c in cList:
        c_path = main_dir+c
        for sub in os.listdir(c_path):
            s = sub.split('_')[0]
            if s in overlapList:
                sub_dir=c_path+'/'+sub
                job_name=s+'_brc_bqsr'
                error=sub_dir+'/'+s+'_brc_bqsr.e'
                out=sub_dir+'/'+s+'_brc_bqsr.o'
                slurm_script=sub_dir+'/'+s+'_brc_bqsr.slurm'
                #tmpdir=sub_dir+'/TmpDir'
                in_bam=sub_dir+'/'+s+'_md_srt.bam'
                out_table=sub_dir+'/'+s+'_recal_data.table'
                out_bam=sub_dir+'/'+s+'_recal.bam'
#                os.system("mkdir -p {0}".format(tmpdir))
                os.system("echo -e '#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name={0}\n#SBATCH --partition=Orion\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=4\n#SBATCH --mem=150gb\n#SBATCH --time=1-00:00:00\n#SBATCH --error={1}\n#SBATCH --output={2}\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=kkingwoo@uncc.edu\n#################################################\n\nmodule load gatk/4.2.0\n\ngatk BaseRecalibrator -I {3} -R {4} --known-sites {5} -O {6}\n\ngatk ApplyBQSR -R {4} -I {3} --bqsr-recal-file {6} -O {7}\n\nrm {3}' > {8}".format(job_name,error,out,in_bam,ref,known,out_table,out_bam,slurm_script))

                os.system("sbatch {0}".format(slurm_script))

overlapped.close()
