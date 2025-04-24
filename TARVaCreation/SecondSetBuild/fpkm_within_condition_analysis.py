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

def parse_stringtie(strip,rdict,sdict,odict,db):
    tpm_dict = {}
    out_dict = {}
    tab_lists = []
    final_results = []
    for condition in sdict.keys():
        tpm_dict[condition] = {}
        with futures.ProcessPoolExecutor(max_workers=32) as mst:
            print(f'************Getting fpkm values for {condition}',datetime.now())
            wait_for = [mst.submit(ReadMstrg.read_string_fpkm,condition,strip,sample,odict) for sample in sdict[condition] ]
            for fu in futures.as_completed(wait_for):
                current = fu.result()
                final_results.append(current)
    for fin in final_results:
    
        for fi in fin:
            tab_lists.append(fi)
    
    dbcon = sql.connect(db)
    df = pd.DataFrame(tab_lists, columns=['ensg', 'enst', 'samp','con','fpkm'])
    df.to_sql("FPKM_PerSample_tab",dbcon,if_exists="replace")

    return 

def build_table(ods,tfile):
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

def get_vep_preds_dtu(vfile,db,tab):
    col_names = ['Gene','Consequence','Feature']
    final_df = pd.DataFrame(columns = col_names) 
    ad_out,control_out = open('ad_vep_cons.csv','w'),open('control_vep_cons.csv','w')
    ad_write,control_write = csv.writer(ad_out),csv.writer(control_out)
    out_dict = {}
    out_dict['AD'] = {}
    out_dict['Control'] = {}
    ad_list,con_list = [],[]
    zn_dict,fpkm_dict = {},{}
    vep,tabs = pd.read_csv(vfile,delimiter='\t'), open(tab,'r')
    ta_read = csv.reader(tabs)
    next(ta_read)
    dbcon = sql.connect(db)
    d = dbcon.cursor()
    all_info_tab,fpkm_tab = "All_Info_TopGenes_tab","FPKM_PerSample_Tab"
    for t in ta_read:
        ensg1,pos1,hgvs1 = t[0],t[1],t[7]
        ensts = vep.loc[vep['#Uploaded_variation'] == hgvs1,['Gene','Consequence','Feature']]
        if not ensts.empty:
            final_df = pd.concat([final_df,ensts],ignore_index = True)
            enst_list = ensts['Feature'].values.tolist()
            que1 = f"SELECT rids FROM {all_info_tab} WHERE poss = ? AND ref_base = ?"
            q1s = d.execute(que1,(pos1,ensg1)).fetchall()
            for q1 in q1s:
                rid = q1[0]
                que2 = f"SELECT ensg,enst,con,fpkm FROM {fpkm_tab} where samp = ?"
                q2s = d.execute(que2,(rid,)).fetchall()
                for q2 in q2s:
                    ensg2,enst2,condition,fpkm = q2[0],q2[1],q2[2],q2[3]
                    ensg = ensg2.split('_')[0]
                    if ensg == ensg1:
                        en2 = enst2.replace('_','.')
                                   
                        if en2 in enst_list:
                            if ensg == "ENSG00000198040":
                                if not en2 in zn_dict:
                                    zn_dict[en2] = {}
                                if not condition in zn_dict[en2]:
                                    zn_dict[en2][condition] = []
                                zn_dict[en2][condition].append(fpkm)
                            if not ensg in fpkm_dict:
                                fpkm_dict[ensg] = {}
                            if not en2 in fpkm_dict[ensg]:
                                fpkm_dict[ensg][en2] = {}
                            if not condition in fpkm_dict[ensg][en2]:
                                fpkm_dict[ensg][en2][condition] = []
                            fpkm_dict[ensg][en2][condition].append(fpkm)

    print(f"ZNF84 dictionary length --> {len(zn_dict)}")
    
    for gene in fpkm_dict:
        for transcript in fpkm_dict[gene]:
            for cond in fpkm_dict[gene][transcript]:
                fpkm_list = fpkm_dict[gene][transcript][cond]
                tot_samps = len(fpkm_list)
                if not tot_samps <= 3:
                    avg_fpkm = sum(fpkm_list)/tot_samps
                    if avg_fpkm > 0.01:
                        cons = final_df.loc[final_df['Feature'] == transcript,['Consequence']].values.tolist()
                        for con in cons:
                            conseq = con[0]
                            if not conseq in out_dict[cond]:
                                out_dict[cond][conseq] = 0
                            out_dict[cond][conseq]+=1
    for key in out_dict:
        for cons,count in out_dict[key].items():
            out_list = [cons,count]
            if key == 'AD':
                ad_list.append(out_list)
            if key == 'Control':
                con_list.append(out_list)
    for a in ad_list:
        ad_write.writerow(a)
    for c in con_list:
        control_write.writerow(c)
                

            
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
    #print(res[0])
    #strings = parse_stringtie(string,res[0],res[1],res[2],dbp)
    #comps = calc_tpms(strings,tfile)
    #hgvss = get_hgvs(tab3,hgvs_out)
    #dtus = analyze_dtu(tab3,tpms,dtu_sig) 
    vep_preds = get_vep_preds_dtu(vep_in,dbp,tab3)
