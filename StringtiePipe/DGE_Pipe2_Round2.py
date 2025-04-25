import os

ad_path, control_path ='/path/to/RNAseq/AD/', '/path/to/RNAseq/Control/'
gff='gencode.v42.annotation.gff3'
outGtf='All_mrgd.gtf'
os.system("rm {0}".format(outGtf))

gList = []

for da in os.listdir(ad_path):
    da_path=ad_path+da
    for ga in os.listdir(da_path):
        ga_path=da_path+'/'+ga
        if ga_path.endswith('_round1.gtf') and not ga_path.endswith('_cov_round1.gtf'):
            gList.append(ga_path)

for dc in os.listdir(control_path):
    dc_path=control_path+dc
    for gc in os.listdir(dc_path):
        gc_path=dc_path+'/'+gc
        if gc_path.endswith('_round1.gtf') and not gc_path.endswith('_cov_round1.gtf'):
            gList.append(gc_path)


gString=''

for gt in gList:
    gString += str(gt)+' '

os.system("stringtie --merge -G {0} -o {1} {2} ".format(gff,outGtf,gString))

