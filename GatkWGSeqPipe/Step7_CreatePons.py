import csv
import os

intervals = '/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/somatic-hg38_CNV_and_centromere_blacklist.hg38liftover.list'
pons = '/nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/pons.args'
ref = '/users/kkingwoo/FunkLab/BWA_Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
dic = '/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.dict'
gr = '/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/somatic-hg38_af-only-gnomad.hg38.vcf'
slurm='/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/createPons.slurm'
job_name='createPons'
out='/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/createPons.o'
err='/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/createPons.e'

out_vcf='/nobackup/kfunk_research/ROSMAP_WGS/TrimmedFASTQFiles/pons.vcf.gz'
variants='/users/kkingwoo/FunkLab/WGS_pipeScripts_9_6_2022/CombinedForPONS.vcf'
os.system("echo -e '#!/bin/sh\n########### Script Settings ####### Katie Jensen #\n#################################################\n#SBATCH --job-name={0}\n#SBATCH --partition=Orion\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=4\n#SBATCH --mem=150gb\n#SBATCH --time=1-00:00:00\n#SBATCH --error={1}\n#SBATCH --output={2}\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=kkingwoo@uncc.edu\n#################################################\n\nmodule load gatk/4.2.0\n\ngatk CreateSomaticPanelOfNormals -R {3} -L {4} --germline-resource {5} -V {6} -O {7}' > {8}".format(job_name,err,out,ref,intervals,gr,variants,out_vcf,slurm))
os.system("sbatch {0}".format(slurm))
