import os
import csv
import subprocess
import json
import pandas as pd
import sqlite3 as sql
import sys
from runqc import RunQC
from datetime import datetime
from concurrent import futures
import tempfile

def get_meta(mono,clin,bb,mono_fq,bubr_fq):
    #mono_dict50, mono_dict76 = {},{}   
    bb_dict= {}
    
    bb_dict['MCI'] = {}
    
    #status_list = [1.0,2.0,4.0]

    #keylist = ['Control','MCI','AD']
    #for ke in keylist:
    #    mono_dict50[ke] = {}
    #    mono_dict76[ke] = {}

    #monocyte_dat = pd.read_csv(mono,delimiter='\t')
    #mono_samps = monocyte_dat['individualID'].unique().tolist()

    bulkbrain_dat = pd.read_csv(bb,delimiter='\t')
    bb_samps = bulkbrain_dat['individualID'].unique().tolist()
    
    info = pd.read_csv(clin)
    
    #monos = info.loc[info['individualID'].isin(mono_samps),["msex","educ","race","spanish","apoe_genotype","age_at_visit_max","age_first_ad_dx","age_death","cts_mmse30_first_ad_dx","cts_mmse30_lv","pmi","braaksc","ceradsc","cogdx","dcfdx_lv","individualID"]]
    #filt_monos = monos.loc[monos['dcfdx_lv'].isin(status_list)]
    
    bbs = info.loc[info['individualID'].isin(bb_samps),["msex","educ","race","spanish","apoe_genotype","age_at_visit_max","age_first_ad_dx","age_death","cts_mmse30_first_ad_dx","cts_mmse30_lv","pmi","braaksc","ceradsc","cogdx","dcfdx_lv","individualID"]]
    filt_bbs = bbs.loc[bbs['dcfdx_lv']==2.0]

    # Lists each value and the number of samples for each 1.0 == Control; 2.0 == MCI; 4.0 == AD
    #mono_con_cts = filt_monos['dcfdx_lv'].value_counts()

    #mfqs,bfqs = [m for m in os.listdir(mono_fq) if os.path.isfile(os.path.join(mono_fq, m)) and m.endswith('fastq.gz')], [b for b in os.listdir(bubr_fq) if os.path.isfile(os.path.join(bubr_fq, b)) and b.endswith('fastq.gz')]
    
    bfqs = [b for b in os.listdir(bubr_fq) if os.path.isfile(os.path.join(bubr_fq, b)) and b.endswith('fastq.gz')]
    
    #for mf in mfqs:
    #    mfs = monocyte_dat.loc[monocyte_dat['name'] == mf, ['individualID','readLength']].values
    #    rid,length = mfs[0][0],mfs[0][1]
    #    mat = filt_monos.loc[filt_monos['individualID'] == rid, 'dcfdx_lv'].values
    #    if mat.size > 0:
    #        con = float(mat[0])
    #        idx = status_list.index(con)
    #        key = keylist[idx]
    #        if length == 50:
    #            if not rid in mono_dict50[key].keys():
    #                mono_dict50[key][rid] = []
    #            fi = mono_fq+mf
    #            mono_dict50[key][rid].append(fi)
    #        if length == 76:
    #            if not rid in mono_dict76[key].keys():
    #                mono_dict76[key][rid] = []
    #            fi = mono_fq+mf
    #            mono_dict76[key][rid].append(fi)
    for bf in bfqs:
        bfs = bulkbrain_dat.loc[bulkbrain_dat['name'] == bf, 'individualID'].values
        mat = filt_bbs.loc[filt_bbs['individualID'] == bfs[0], 'dcfdx_lv'].values
        if mat.size > 0:
            rid = bfs[0]
            if not rid in bb_dict['MCI'].keys():
                bb_dict['MCI'][rid] = []
            fil = bubr_fq+bf
            bb_dict['MCI'][rid].append(fil)
    
    return bb_dict  #mono_dict50,mono_dict76,bb_dict

#def trim(md50,md76):  #,bd150):
def trim(bd150):
    final_results = []
    #specs_list = ('50','76') #,'150')
    specs_list = '150'
    #tissues_list = ('monocytes','monocytes') #,'bulkbrain')
    tissues_list = 'bulkbrain'
    #dict_list = (md50,md76) #,bd150)
    dict_list = bd150
    
    adapter_file = 'trimmomatic-0.39-2/adapters/TruSeq3-PE.fa'
    with futures.ProcessPoolExecutor(max_workers=8) as mst:
        print('************ Running QC on the files --> ',datetime.now())
        wait_for = [mst.submit(RunQC.run_qc,dict_list,specs_list,tissues_list,adapter_file)]
        #wait_for = [mst.submit(RunQC.run_qc,dict_list[n],specs_list[n],tissues_list[n],adapter_file) for n in range(0,len(dict_list))]
        for fu in futures.as_completed(wait_for):
            current = fu.result()
            final_results.append(current)
    out_dict = {}
    #out_dict['monocytes'] ={}
    out_dict['bulkbrain'] = {}
    #out_dict['monocytes']['50'] = []
    #out_dict['monocytes']['76'] = []
    out_dict['bulkbrain']['150'] = []
    for fin in final_results:
        for fi in fin:
            if len(fi) > 1:
                read_len = fi[0]
                for tupe in range(1,len(fi)):
                    tu = fi[tupe]
                    #if read_len == 50:
                    #    out_dict['monocytes']['50'].append(tu)
                    #if read_len == 76:
                    #    out_dict['monocytes']['76'].append(tu)
                    if read_len == 150:
                        out_dict['bulkbrain']['150'].append(tu)
    
    return out_dict

def star(diction,gtf,gen,sl,dfile,gatk,dout):
    with open(dfile,'w') as df:
        dumped = json.dump(diction,df)
        df.close()
        fa = "GRCh38.primary_assembly.genome.fa"
        ##Change the following in the 'slurm' variable before running STAR: mem=175gb; cpus-per-task=16; time=4-00:00:00
        slurm = f"""#!/bin/sh
########### Script Settings ####### Katie Jensen #
#################################################
#SBATCH --job-name=mm_star
#SBATCH --partition=
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --error="star_par.e"
#SBATCH --output="star_par.o"
#################################################

module load star

GTF="{gtf}"
GENOME="{gen}"
DICT="{dfile}"
FASTA="{fa}"
GATK_SLURM="{gatk}"
DUMPOUT="{dout}"

python star.py "$GENOME" "$GTF" "$DICT" "$SLURM_CPUS_PER_TASK" "$GATK_SLURM" "$DUMPOUT"

"""
    
    with open(sl,'w') as file:
        file.write(slurm)

        file.close()
    print(f"SLURM script written to {sl}")
    run_job = f"sbatch {sl}"
    result = subprocess.run(run_job, shell=True, check=True)
    print(result.stdout)
    print(result.stderr)
    return

if __name__=='__main__':
    dbp = sys.argv[1]
    mon_dat = sys.argv[2]
    clin = sys.argv[3]
    bb_dat = sys.argv[4]
    bb_fast = sys.argv[5]
    mon_fast = sys.argv[6]
    gt=sys.argv[7]
    geno=sys.argv[8]
    slu = sys.argv[9]
    dufi = sys.argv[10]
    g_slurm = sys.argv[11]
    dump_out = sys.argv[12]
    res = get_meta(mon_dat,clin,bb_dat,mon_fast,bb_fast)
    #qc = trim(res[0],res[1],res[2])
    qc = trim(res)
    star_submit = star(qc,gt,geno,slu,dufi,g_slurm,dump_out)
