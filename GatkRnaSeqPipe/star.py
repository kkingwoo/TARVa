import sys
import json
import subprocess
from datetime import datetime
from concurrent import futures
from star_run import RunSTAR


def align(gen,gt,dic,cpus_t,dumpout):
    downstream_sets = []
    final_results = []
    di = None
    with open(dic, 'r') as dumpfile:
        di = json.load(dumpfile)
    tissue_list,read_len_list,tupe_list = [],[],[]
    for tissue in di.keys():
        for readlen in di[tissue].keys():
            li = di[tissue][readlen]
            for tup in li:
                read_len_list.append(readlen)
                tissue_list.append(tissue)
                tupe_list.append(tup)
    with futures.ProcessPoolExecutor(max_workers=4) as mst:
        print('************ Starting STAR in parallel --> ',datetime.now())
        wait_for = [mst.submit(RunSTAR.run_star,tupe_list[i],gen,tissue_list[i],read_len_list[i],cpus_t) for i in range(0,len(tupe_list))]
        for fu in futures.as_completed(wait_for):
            current = fu.result()
            final_results.append(current)
    for fi in final_results:
        if fi:
            downstream_sets.append(fi)
    with open(dumpout,'w') as jout:
        json.dump(downstream_sets,jout)
    return dumpout

def gatk_slurm(tupes,slurm,dumpout):
    ds_var = "monothreecon_bulkbrainmci"
    fa_string = "GRCh38.primary_assembly.genome.fa"   
    known_file = "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
    sl = f"""#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=full_gatk
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=128gb
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --error="full_gatk.e"
#SBATCH --output="full_gatk.o"
#################################################

module load anaconda3
source activate bio_qc
module load R/4.3.3
export R_LIBS_USER=Rlibs

DUMPED="{dumpout}"
FA="{fa_string}"
DS_VAR="{ds_var}"
KNOWN="{known_file}"

python full_gatk.py "$DUMPED" "$FA" "$DS_VAR" "$KNOWN"

"""
    with open (slurm,'w') as outslurm:
        outslurm.write(sl)
        outslurm.close()

    print(f"SLURM script written to {slurm}")
    run_job = f"sbatch {slurm}"
    result = subprocess.run(run_job, shell=True, check=True)
    print(result.stdout)
    print(result.stderr)
    return



if __name__=='__main__':
    genome = sys.argv[1]
    gtf = sys.argv[2]
    dumpy = sys.argv[3]
    cpus = int(sys.argv[4])
    slu = sys.argv[5]
    dumpo = sys.argv[6]
    aln = align(genome,gtf,dumpy,cpus,dumpo)
    gatk = gatk_slurm(aln,slu,dumpo)



