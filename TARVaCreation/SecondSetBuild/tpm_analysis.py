import os
import sys
import sqlite3 as sql
import pandas as pd
import numpy as np
from datetime import datetime
from concurrent import futures
from read_mstrg_tab import ReadMstrg
from scipy import stats
import csv
import math

def retrieve_tids(db):
    ref_dict,samps_dict,other_dict = {},{},{}
    con = sql.connect(db)
    c = con.cursor()
    all_tab,samp_tab,fast_tab = 'All_Info_TopGenes_tab','sample_tab','fasta_tab'
    first = c.execute(f"SELECT DISTINCT refid,poss FROM {all_tab}").fetchall()
    for fi in first:
        ref_id,pos = fi[0],fi[1]
        if not ref_id in ref_dict.keys():
            ref_dict[ref_id] = []
        if not pos in ref_dict[ref_id]:
            ref_dict[ref_id].append(pos)
    second = c.execute(f"SELECT DISTINCT condition,rid FROM {samp_tab}")
    for se in second:
        cond,rids = se[0],se[1]
        if not cond in samps_dict.keys():
            samps_dict[cond] = []
        if not rids in samps_dict[cond]:
            samps_dict[cond].append(rids)
    for refs in ref_dict.keys():
        for position in ref_dict[refs]:
            third = f"SELECT DISTINCT ref_id,tid FROM {fast_tab} WHERE ref_id = ? AND pos = ?"
            thir = c.execute(third, (refs,position)).fetchall()
            for thi in thir:
                refids,tids = thi[0],thi[1]
                if not refids in other_dict.keys():
                    other_dict[refids] = []
                if not tids in other_dict[refids]:
                    other_dict[refids].append(tids)
    return ref_dict,samps_dict,other_dict

def parse_stringtie(strip,rdict,sdict,odict):
    tpm_dict = {}
    out_dict = {}
    final_results = []
    for condition in sdict.keys():
        tpm_dict[condition] = {}
        with futures.ProcessPoolExecutor(max_workers=32) as mst:
            print(f'************Getting tpm values for {condition}',datetime.now())
            wait_for = [mst.submit(ReadMstrg.read_string_tpm,condition,strip,sample,odict) for sample in sdict[condition] ]
            for fu in futures.as_completed(wait_for):
                current = fu.result()
                final_results.append(current)
    for fin in final_results:
        for tu in fin:
            condi,refs,tis,lens,vals = tu[0],tu[1],tu[2],tu[3],tu[4]
            if not refs in out_dict.keys():
                out_dict[refs] = {}
            if not tis in out_dict[refs].keys():
                out_dict[refs][tis] = {}
            if not condi in out_dict[refs][tis].keys():
                out_dict[refs][tis][condi] = []
            out_dict[refs][tis][condi].append(vals)
    
    return out_dict

def calc_tpms(ods,tfile):
    outfile = open(tfile,'w')
    out_write = csv.writer(outfile)
    header = ['refid','tid','ustat','pval','ad_tot_tpm','con_tot_tpm']
    out_write.writerow(header)
    for key in ods.keys():
        for kkey in ods[key].keys():
            ad,con = [],[]
            out_row = [key]
            out_row.append(kkey)
            for conds in ods[key][kkey].keys():
                if conds == 'AD':
                    ad = ods[key][kkey][conds]
                if conds == 'Control':
                    con = ods[key][kkey][conds]
            ads,cons = sum(ad),sum(con)  
            stat,p_value = stats.mannwhitneyu(ad,con)
            #if float(p_value) < float(0.05):
            out_row.append(stat)
            out_row.append(p_value)
            out_row.append(ads)
            out_row.append(cons)
            out_write.writerow(out_row)

    return

def get_hgvs(t3,hout):
    ho = open(hout,'w')
    tab3 = pd.read_csv(t3)
    ids = tab3['hgvs'].tolist()

    for i in ids:
        ho.write(f'{i}\n')

    return

def analyze_dtu(t3,tpm,dtu):
    
    t3s,ts = pd.read_csv(t3),pd.read_csv(tpm)
    
    refids = t3s['ensg'].unique().tolist()
    sigs = ts[(ts['pval'] < 0.05) & (ts['refid'].str.contains('|'.join(refids) + '_'))]
    
    sigs.to_csv(dtu,index=False)

    return dtu

def get_vep_preds_dtu(vfile,dtu):
    vep = pd.read_csv(vfile,delimiter='\t')

    return

if __name__=='__main__':
    dbp = sys.argv [1]
    string = sys.argv[2]
    tab3 = sys.argv[3]
    tpms = sys.argv[4]
    hgvs_out = sys.argv[5]
    vep_in = sys.argv[6]
    dtu_sig = sys.argv[7]
    #res = retrieve_tids(dbp)
    #strings = parse_stringtie(string,res[0],res[1],res[2])
    #comps = calc_tpms(strings,tfile)
    #hgvss = get_hgvs(tab3,hgvs_out)
    dtus = analyze_dtu(tab3,tpms,dtu_sig) 
    vep_preds = get_vep_preds_dtu(vep_in,dtu_sig)
