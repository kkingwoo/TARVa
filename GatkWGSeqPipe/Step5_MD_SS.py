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
                job_name=s+'_md_ss'
                error=sub_dir+'/'+s+'_md_ss.e'
                out=sub_dir+'/'+s+'_md_ss.o'
                slurm_script=sub_dir+'/'+s+'_md_ss.slurm'
                tmpdir=sub_dir+'/TmpDir'
                ub=sub_dir+'/'+s+'_marked_srtd.bam'
                sb=sub_dir+'/'+s+'_md_srt.bam'
                ab=sub_dir+'/'+s+'_SortedBams/'+s+'_marked_dupes.bam'
                ob=sub_dir+'/'+s+'_SortedBams/'+s+'_mrgd_aln.bam'
                mt=sub_dir+'/'+s+'_dupes.txt'              
                os.system("mkdir -p {0}".format(tmpdir))
                os.system("echo -e '#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name={0}\n#SBATCH --partition=Orion\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=16\n#SBATCH --mem=64gb\n#SBATCH --time=1-00:00:00\n#SBATCH --error={1}\n#SBATCH --output={2}\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=kkingwoo@uncc.edu\n#################################################\n\nmodule load gatk/4.2.0\n\nexport _JAVA_OPTIONS=-Djava.io.tmpdir={3}\n\ngatk MarkDuplicates -I {4} -O {5} -M {6}\n\nrm {4}\n\ngatk SortSam -I {5} -O {7} -R {8} -SO coordinate --CREATE_INDEX true --TMP_DIR {3}\n\nrm {5}' > {9}".format(job_name,error,out,tmpdir,ob,ab,mt,sb,ref,slurm_script))

                os.system("sbatch {0}".format(slurm_script))

overlapped.close()
