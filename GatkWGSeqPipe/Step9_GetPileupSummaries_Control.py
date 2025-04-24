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
gr = '/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/somatic-hg38_small_exac_common_3.hg38.vcf.gz'
pons='/nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/pons.vcf.gz'

main_dir = '/nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/Control/'

for a in os.listdir(main_dir):
    ai = a.split('_')[0]
    if ai in overlapList:
        a_path = main_dir+a
        slurm=a_path+'/'+ai+'_PileupSum.slurm'
        job_name=ai+'_PileupSum'
        out=a_path+'/'+ai+'_PileupSum.o'
        err=a_path+'/'+ai+'_PileupSum.e'
        in_bam=a_path+'/'+ai+'_recal.bam'
        out_table=a_path+'/'+ai+'_pileupSum.table'
        cont_table=a_path+'/'+ai+'_contamination.table'
        out_vcf=a_path+'/'+ai+'_contaminates_filtered.vcf.gz'
        in_vcf=a_path+'/'+ai+'_pon.vcf.gz'

        os.system("echo -e '#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name={0}\n#SBATCH --partition=Orion\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=4\n#SBATCH --mem=150gb\n#SBATCH --time=1-00:00:00\n#SBATCH --error={1}\n#SBATCH --output={2}\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=kkingwoo@uncc.edu\n#################################################\n\ngatk GetPileupSummaries -I {3} -V {4} -L {4} -O {6}\n\ngatk CalculateContamination -I {6} -O {7}\n\ngatk FilterMutectCalls -R {8} -V {5} --contamination-table {7} -O {9}' > {10}".format(job_name,err,out,in_bam,gr,in_vcf,out_table,cont_table,ref,out_vcf,slurm))
        #os.system("module load anaconda3")
        #os.system("conda activate gatk4")
        os.system("sbatch {0}".format(slurm))
