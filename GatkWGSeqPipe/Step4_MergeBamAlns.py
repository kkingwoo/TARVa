import os 
import csv

ref = '/users/kkingwoo/FunkLab/BWA_Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
#ref_dict=${ref%%.fa}".dict"

#gatk CreateSequenceDictionary -R $ref -O $ref_dict

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
                job_name=s
                error=sub_dir+'/'+s+'_mba.e'
                out=sub_dir+'/'+s+'_mba.o'
                slurm_script=sub_dir+'/'+s+'_mba.slurm'
                tmpdir=sub_dir+'/TmpDir'
                ub=sub_dir+'/'+s+'_unaln_mrgd_srtd.bam'
                ab=sub_dir+'/'+s+'_SortedBams/'+s+'_merged.bam'
                ob=sub_dir+'/'+s+'_SortedBams/'+s+'_mrgd_aln.bam'
                
                os.system("mkdir -p {0}".format(tmpdir))
                os.system("echo -e '#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name={0}\n#SBATCH --partition=Orion\n#SBATCH --nodes=1\n#SBATCH --mem=150gb\n#SBATCH --time=1-00:00:00\n#SBATCH --error={1}\n#SBATCH --output={2}\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=kkingwoo@uncc.edu\n#################################################\n\nmodule load gatk/4.2.0\n\nexport _JAVA_OPTIONS=-Djava.io.tmpdir={3}\n\ngatk MergeBamAlignment --TMP_DIR {3} -ALIGNED {4} -UNMAPPED {5} -O {6} -R {7} \n\nrm {5}\n\nrm {4}' > {8}".format(job_name,error,out,tmpdir,ab,ub,ob,ref,slurm_script))
                os.system("sbatch {0}".format(slurm_script))

overlapped.close()
