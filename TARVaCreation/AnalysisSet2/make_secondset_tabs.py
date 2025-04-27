import os
import sqlite3 as sql
from datetime import datetime
import sys
import pandas as pd
import numpy as np

class SampleTabs:
    @staticmethod
    def define_condition(cog_val):
        val_list, return_list = [1.0,2.0,4.0],['Control','MCI','AD']
        idx = val_list.index(cog_val)
        con = return_list[idx]
        return con


    @staticmethod
    def sample_tabs(pa,fil,tissue,clinical,dbc):
        fa_tab = 'fasta_tab'
        fi = os.path.join(pa,fil)
        samp = fi.split('/')[-1].split('_')[0]
        head_list = ["chrom","pos","id","ref","alt","qual","filter","info","format",samp]
        clin_info = pd.read_csv(clinical)
        samp_dat = clin_info.loc[clin_info['individualID'] == samp,["msex","age_death","apoe_genotype","braaksc","ceradsc","dcfdx_lv"]].iloc[0].tolist()
        
        cons = samp_dat[-1]
        condition = SampleTabs.define_condition(cons)
        
        sex,age,apoe,braak,cerad = samp_dat[0], samp_dat[1], samp_dat[2], samp_dat[3], samp_dat[4],
        rf = pd.read_csv(fi,sep='\t',comment='#',header=0,names=head_list)
        rf['condition'],rf['tissue'],rf['rid'],rf['sex'],rf['age'],rf['apoe'],rf['braak'],rf['cerad'],rf['string_id'],rf['ref_id'],rf['tid'],rf['exon_pos'],rf['trans_pos'],rf['a1_edit_type'],rf['a2_edit_type'],rf['a1_prop_exon'],rf['a2_prop_exon'],rf['a1_prop_trans'],rf['a2_prop_trans'],rf['alt1'],rf['alt2'],rf['exon'],rf['a1_ad'],rf['a2_ad'],rf['dp'],rf['per_exon_pos_wt'],rf['per_trans_pos_wt'],rf['exon_len'],rf['trans_len'],rf['strand'] = condition,tissue,samp,sex,age,apoe,braak,cerad,None,None,None,None,None,None,None,np.nan,np.nan,np.nan,np.nan,None,None,None,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,None
        index_names = rf[ ~(rf["chrom"].str.startswith('chr'))].index
        rf.drop(index_names,inplace=True)
        indexes = rf[ (rf['chrom'].str.endswith('Y'))].index
        rf.drop(indexes,inplace=True)
        nuc_set = {'A','T','C','G'}
        
        db_con = sql.connect(dbc,check_same_thread=False)
        db = db_con.cursor()
        out_list = []
        for index,row in rf.iterrows():
            chrom = row['chrom']
            pos = row['pos']
            ref = row['ref']
            read_counts = row[samp]
            alts = row['alt']
            rc = read_counts.split(':')[1:3]
            dp = int(rc[1])
           ##Filter for depth of reads at each position ##
            if dp >= 20:
                re = ref
                ref_one = re[0]
                al = alts.split(',')
                alt1 = al[0]
               ## Filter out positions where there is no alternate nucleotide call ##
                if alt1 != '.':
                    a1_prop_exon,a1_prop_trans,a2_prop_exon,a2_prop_trans = float(),float(),float(),float()
                    type1 = None
                    alt2 = None
                    type2 = None
                    ads = rc[0].split(',')
                    ref_ad = int(ads[0])
                    a1_ad,a2_ad = int(),int()
                    single_iso = ''
                    ## Where there is only one alternate nucleotide ##
                    if not ',' in alts:
                        a1_ad = int(ads[1])
                    ## Where there are two alternate nucleotides ##
                    if ',' in alts:
                        alt2 = al[1]
                        a1_ad,a2_ad = int(ads[1]),int(ads[2])
                    gets = f"SELECT * FROM {fa_tab} WHERE chrom = ? AND pos = ? and nuc = ?"
                    get = db.execute(gets,(chrom,pos,ref_one)).fetchall()
                    for g in get:
                        strand,string_id,ref_id,tid,exon_pos,trans_pos,pep,ptp,exon,exon_len,trans_len,per_exon_pos_wt,per_trans_pos_wt = g[4],g[12],g[13],g[2],g[5],g[6],g[9],g[10],g[11],g[7],g[8],float(g[9]),float(g[10])
                        a1_prop_exon,a1_prop_trans = float((a1_ad/dp) * pep),float((a1_ad/dp) * ptp)
                        if alt2 != None:
                            a2_prop_exon,a2_prop_trans = float((a2_ad/dp) * pep),float((a2_ad/dp) * ptp)
       
                        rf.at[index,'string_id'] = string_id
                        rf.at[index,'ref_id'] = ref_id.split('_')[0]
                        rf.at[index,'tid'] = tid
                        rf.at[index,'exon_pos'] = exon_pos
                        rf.at[index,'trans_pos'] = trans_pos
                        rf.at[index,'a1_edit_type'] = type1
                        rf.at[index,'a2_edit_type'] = type2
                        rf.at[index,'a1_prop_exon'] = a1_prop_exon
                        rf.at[index,'a2_prop_exon'] = a2_prop_exon
                        rf.at[index,'a1_prop_trans'] = a1_prop_trans
                        rf.at[index,'a2_prop_trans'] = a2_prop_trans
                        rf.at[index,'alt1'] = alt1
                        rf.at[index,'alt2'] = alt2
                        rf.at[index,'exon'] = exon
                        rf.at[index,'a1_ad'] = a1_ad
                        rf.at[index,'a2_ad'] = a2_ad
                        rf.at[index,'dp'] = dp
                        rf.at[index,'strand'] = strand
                        rf.at[index,'exon_len'] = exon_len
                        rf.at[index,'trans_len'] = trans_len
                        rf.at[index,'per_exon_pos_wt'] = per_exon_pos_wt
                        rf.at[index,'per_trans_pos_wt'] = per_trans_pos_wt

                        
        drop_cols = ['id','alt','qual','filter','format','info',samp]
        rf.drop(columns=drop_cols,inplace=True)
        
        return rf

    

