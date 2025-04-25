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

slurm='importDB_forPons.slurm'
job_name='importDB_forPons'
out='importDB_forPons.o'
err='importDB_forPons.e'
os.system("echo -e '#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name={0}\n#SBATCH --partition=\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=4\n#SBATCH --mem=150gb\n#SBATCH --time=1-00:00:00\n#SBATCH --error={1}\n#SBATCH --output={2}\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=\n#################################################\n\nmodule load gatk/4.2.0\n\ngatk GenomicsDBImport -R {3} -L {4} --genomicsdb-workspace-path pon_db -V {5}' > {6}".format(job_name,err,out,ref,intervals,pons,slurm))
os.system("sbatch {0}".format(slurm))
