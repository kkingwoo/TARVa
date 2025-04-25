import os

os.system("module load gatk/4.2.0")
intervals = 'somatic-hg38_CNV_and_centromere_blacklist.hg38liftover.list'
pons = '/path/to/ROSMAP_WGS/TrimmedFASTQFiles/pons.args'
ref = 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
#pList = []
#for p in pons.readlines():
#    pList.append(p.strip('\n'))

#string1 = pList[0]
#string = '-V '.join(pList)

slurm='SelectVars_forPons.slurm'
job_name='SelectVars_forPons'
out='SelectVars_forPons.o'
err='SelectVars_forPons.e'
pons_vcf='CombinedForPONS.vcf'
os.system("echo -e '#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name={0}\n#SBATCH --partition=\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=4\n#SBATCH --mem=150gb\n#SBATCH --time=1-00:00:00\n#SBATCH --error={1}\n#SBATCH --output={2}\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=\n#################################################\n\nmodule load gatk/4.2.0\n\ngatk SelectVariants -R {3} -V gendb:///nobackup/kfunk_research/ROSMAP_WGS/pon_db -O {4}' > {5}".format(job_name,err,out,ref,pons_vcf,slurm))
os.system("sbatch {0}".format(slurm))

