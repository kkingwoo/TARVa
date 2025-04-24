import os
import csv
import shutil

main_dir = '/projects/kfunk_research/ROSMAP/RNAseq/ROSMAP_BulkBrain/ROSMAP_BulkBrain_Trimmed_fastq/'

out_file = main_dir+'ROSMAP_RNAseq_designMat.txt'
design_mat = open(out_file, 'w')
mat_write = csv.writer(design_mat, delimiter='\t')

header = ["ID", "condition"]
mat_write.writerow(header)

cList = ["AD", "Control"]

R_Dir = main_dir+"DGEFiles_ForStringtie"
os.mkdir(R_Dir)

for d in os.listdir(main_dir):
    if d in cList:
        d_path = main_dir+d
        for subdir in os.listdir(d_path):
            r_sub = R_Dir+'/'+subdir
            os.mkdir(r_sub)
            w = [subdir, d]
            mat_write.writerow(w)
            f_path = d_path+'/'+subdir
            tList = []
            for t in os.listdir(f_path):
                if t.endswith('ctab'):
                    tList.append(t)
            if len(tList) == 0:
                slurm = f_path+'/'+subdir+'_dge2_3.slurm'
                os.system("sbatch {0}".format(slurm))
            else:
                for tab in tList:
                    old = f_path+'/'+tab
                    new = r_sub+'/'+tab
                    print(old, new)
                    shutil.move(old, new)
