import subprocess
import os
import re
import csv
import json

class RunSTAR:
    @staticmethod
    def run_star(fss,genome,tiss,readlen,cpus_task):
        f1,f2 = '',''
        f_pref = fss[0]
        fsplit = f_pref.split('/')
        outpath = '/'.join(fsplit[:-2])+'/'
        rid = fsplit[-1].split('_')[0]
        permissions = 0o770
        out_pref = f"{outpath}{rid}_{tiss}_{readlen}"
        out1 = os.path.join(out_pref,"2PassMode/")
        os.makedirs(out1,exist_ok=True)
        os.chmod(out1,permissions)
        
        #tab = os.path.join(out1,"SJ.out.tab")
        for fs in fss:
            if "END1" in fs or "R1_001" in fs:
                f1 = fs
            if "END2" in fs or "R2_001" in fs:
                f2 = fs

        #string_state = f"STAR --runThreadN {cpus_task} --genomeDir {genome} --readFilesCommand zcat  --readFilesIn {f1} {f2} --outFileNamePrefix {out1} --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outSAMmode Full --outSAMattributes NH HI AS nM XS --outSAMstrandField intronMotif --limitBAMsortRAM 30000000000 " 
        #subprocess.run(string_state, shell=True,check=True)
        
        #outs = out1+"Aligned.sortedByCoord.out.bam"
        bam = f"{out1}{rid}_{tiss}_{readlen}_2Pass.bam"
        #ob = os.system(f"mv {outs} {bam}")
        log_file = f"{out1}Log.final.out"
        out_tupe = ()
        li = []
    
        logs = open(log_file,'r')
        for log in logs.readlines():
            lo = log.split('\t')
            for l in lo:
                if '                        Uniquely mapped reads % |' == l:
                    perc = float(lo[-1].split('\n')[0].split('%')[0])
                    if perc >= 75.00:
                        out_tupe = (f1,f2,bam,rid)
                
        
        return out_tupe


