import os
import vcf
import csv


rna_path ='/path/to/rnaseq/vcfs/'
wgs_path ='/path/to/wgs/vcfs/' 

wDirDict = {}

for di in os.listdir(wgs_path):
    w_con_path = wgs_path+di+'/'
    if os.path.isdir(w_con_path):
        for subdir in os.listdir(w_con_path):
            s = subdir.split('_')[0]
            wDirDict[s] = ''
            subdir_path = w_con_path+subdir+'/'
            for fi in os.listdir(subdir_path):
                if fi.endswith("contaminates_filtered.vcf.gz"):
                    fi_path = subdir_path+fi
                    wDirDict[s] = fi_path
matchedDirDict = {}

for d in os.listdir(rna_path):
    r_con_path = rna_path+d+'/'
    if os.path.isdir(r_con_path):
        for sd in os.listdir(r_con_path):
            if sd in wDirDict.keys():
                r_samp_path = r_con_path+sd+'/'
                for f in os.listdir(r_samp_path):
                    if f.endswith("varcall_flt.vcf"):
                        f_path = r_samp_path+f
                        w_path = wDirDict[sd]
                        matchedDirDict[sd] = []
                        matchedDirDict[sd].append(f_path)
                        matchedDirDict[sd].append(w_path)
                    

for k in matchedDirDict.keys():
    out_dir = '/'.join(matchedDirDict[k][0].split('/')[0:-1])+'/'
    temp_varDict = {}
    temp_varDict['WGS'] = {}
    rna_vcf = matchedDirDict[k][0] 
    wgs_vcf = matchedDirDict[k][1]
    os.system('gunzip {0}'.format(wgs_vcf))
    wgs_vcf = wgs_vcf.split('.')[0]+'.vcf'
    wgs = vcf.Reader(filename = wgs_vcf)
    rna = vcf.Reader(filename = rna_vcf)
    nombre = out_dir+k+'_rna_unique.vcf'
    rna_unique = vcf.Writer(open(nombre, 'w'), rna)
    for record in wgs:
        chrom, pos, ref, alt1 = record.CHROM, record.POS, record.REF, record.ALT
        varcall = str(ref)+'-'+str(alt1[0])
            
        if not chrom in temp_varDict['WGS'].keys():
            temp_varDict['WGS'][chrom] = {}
            temp_varDict['WGS'][chrom][pos] = ''
            temp_varDict['WGS'][chrom][pos] = varcall
        else:
            if not pos in temp_varDict['WGS'][chrom].keys():
                temp_varDict['WGS'][chrom][pos] = ''
                temp_varDict['WGS'][chrom][pos] = varcall
            else:
                temp_varDict['WGS'][chrom][pos] = varcall
    for rec in rna:
        chro, po, re, alt = rec.CHROM, rec.POS, rec.REF, rec.ALT 
        vc = str(re)+'-'+str(alt[0])
        if chro in temp_varDict['WGS'].keys():
            if not po in temp_varDict['WGS'][chro].keys():
                rna_unique.write_record(rec)
                    

