import os
import sqlite3 as sql
from datetime import datetime
import sys
import pandas as pd
import numpy as np

def parse_dict_list(diction,pos):
    for key, value in diction.items():
        if isinstance(value,set) and pos in value:
            return key
    return None

def a1_pos_snp(rf,a1,t1):
    edit = rf+'-'+a1
    if edit == 'A-G':
        t1 = 'A-I'
    elif edit == 'C-T':
        t1 = 'C-U'
    elif a1 == '*':
        t1 = 'del'
    else:
        t1 = 'other'
    return t1

def a1_neg_snp(rf,a1,t1):
    edit = rf+'-'+a1
    if edit == 'T-C':
        t1 = 'A-I'
    elif edit == 'G-A':
        t1 = 'C-U'
    elif a1 == '*':
        t1 = 'del'
    else:
        t1 = 'other'
    return t1

def a1_pos_indel(rf,a1,t1):
    ref_one = rf[0]
    a1_one = a1[0]
    edit = ref_one+'-'+a1_one
    if len(rf) < len(a1):
        if edit == 'A-G':
            t1 = 'A-I_insert'
        elif edit == 'C-T':
            t1 = 'C-U_insert'
        else:
            t1 = 'insert'
    else:
        if len(rf) > len(a1):
            if edit == 'A-G':
                t1 = 'A-I_del'
            elif edit == 'C-T':
                t1 = 'C-U_del'
            else:
                t1 = 'del'
    return t1

def a1_neg_indel(rf,a1,t1):
    ref_one = rf[0]
    a1_one = a1[0]
    edit = ref_one+'-'+a1_one
    if len(rf) < len(a1):
        if edit == 'T-C':
            t1 = 'A-I_insert'
        elif edit == 'G-A':
            t1 = 'C-U_insert'
        else:
            t1 = 'insert'
    else:
        if len(rf) > len(a1):
            if edit == 'T-C':
                t1 = 'A-I_del'
            elif edit == 'G-A':
                t1 = 'C-U_del'
            else:
                t1 = 'del'
    return t1

def a2_pos_snp(rf,a2,t2):
    edit = rf+'-'+a2
    if edit == 'A-G':
        t2 = 'A-I'
    elif edit == 'C-T':
        t2 = 'C-U'
    elif a2 == '*':
        t2 = 'del'
    else:
        t2 = 'other'
    return t2

def a2_neg_snp(rf,a2,t2):
    edit = rf+'-'+a2
    if edit == 'T-C':
        t2 = 'A-I'
    elif edit == 'G-A':
        t2 = 'C-U'
    elif a2 == '*':
        t2 = 'del'
    else:
        t2 = 'other'
    return t2

def a2_pos_indel(rf,a2,t2):
    ref_one = rf[0]
    a2_one = a2[0]
    edit = ref_one+'-'+a2_one
    if len(rf) < len(a2):
        if edit == 'A-G':
            t2 = 'A-I_insert'
        elif edit == 'C-T':
            t2 = 'C-U_insert'
        else:
            t2 = 'insert'
    else:
        if len(rf) > len(a2):
            if edit == 'A-G':
                t2 = 'A-I_del'
            elif edit == 'C-T':
                t2 = 'C-U_del'
            else:
                t2 = 'del'
    return t2

def a2_neg_indel(rf,a2,t2):
    ref_one = rf[0]
    a2_one = a2[0]
    edit = ref_one+'-'+a2_one
    if len(rf) < len(a2):
        if edit == 'T-C':
            t2 = 'A-I_insert'
        elif edit == 'G-A':
            t2 = 'C-U_insert'
        else:
            t2 = 'insert'
    else:
        if len(rf) > len(a2):
            if edit == 'T-C':
                t2 = 'A-I_del'
            elif edit == 'G-A':
                t2 = 'C-U_del'
            else:
                t2 = 'del'
    return t2


class SampleTabs:
    def instance_parsing(self,diction,pos):
        trans_tab = parse_dict_list(diction,pos)
        return trans_tab

    def assign_func_one(self,rf,a1,t1):
        types = a1_pos_snp(rf,a1,t1)
        return types
    def assign_func_two(self,rf,a1,t1):
        types = a1_pos_indel(rf,a1,t1)
        return types
    def assign_func_three(self,rf,a1,t1):
        types = a1_neg_snp(rf,a1,t1)
        return types
    def assign_func_four(self,rf,a1,t1):
        types = a1_neg_indel(rf,a1,t1)
        return types
    def assign_func_five(self,rf,a2,t2):
        types = a2_pos_snp(rf,a2,t2)
        return types
    def assign_func_six(self,rf,a2,t2):
        types = a2_pos_indel(rf,a2,t2)
        return types
    def assign_func_seven(self,rf,a2,t2):
        types = a2_neg_snp(rf,a2,t2)
        return types
    def assign_func_eight(self,rf,a2,t2):
        types = a2_neg_indel(rf,a2,t2)
        return types

    @staticmethod
    def sample_tabs(condition,samp,matched,dbc,gd):
        fa_tab = 'fasta_tab'
        head_list = ["chrom","pos","id","ref","alt","qual","filter","info","format",samp]
        rf,gf = pd.read_csv(matched[condition][samp][0],sep='\t',comment='#',header=0,names=head_list),pd.read_csv(matched[condition][samp][1],sep='\t',comment='#',header=0,names=head_list)
        rf['condition'],rf['rid'],rf['string_id'],rf['ref_id'],rf['tid'],rf['exon_pos'],rf['trans_pos'],rf['single_iso'],rf['a1_edit_type'],rf['a2_edit_type'],rf['a1_prop_exon'],rf['a2_prop_exon'],rf['a1_prop_trans'],rf['a2_prop_trans'],rf['alt1'],rf['alt2'],rf['exon'],rf['a1_ad'],rf['a2_ad'],rf['dp'],rf['per_exon_pos_wt'],rf['per_trans_pos_wt'],rf['exon_len'],rf['trans_len'],rf['strand'] = condition,samp,None,None,None,None,None,None,None,None,np.nan,np.nan,np.nan,np.nan,None,None,None,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,None
        merged_df = pd.merge(rf,gf,on=['chrom','pos'], how='left',suffixes=('_rna','_genomic'))
        rf['genomic'] = np.where(pd.notna(merged_df['alt_genomic']), 'yes', 'no')
        index_names = rf[ ~(rf["chrom"].str.startswith('chr'))].index
        rf.drop(index_names,inplace=True)
        indexes = rf[ (rf['chrom'].str.endswith('Y'))].index
        rf.drop(indexes,inplace=True)
        nuc_set = {'A','T','C','G'}
        rna_frame = rf[rf['genomic'] == 'no']
        db_con = sql.connect(dbc,check_same_thread=False)
        db = db_con.cursor()
        out_list = []
        for index,row in rna_frame.iterrows():
            chrom = row['chrom']
            pos = row['pos']
            ref = row['ref']
            read_counts = row[samp]
            alts = row['alt']
            rc = read_counts.split(':')[1:3]
            dp = int(rc[1])
            ## Filter for depth of reads at each position ##
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
                        ## Single isoform gene? ##
                        if gd[chrom][string_id] == 1:
                            single_iso = 'yes'
                        if gd[chrom][string_id] > 1:
                            single_iso = 'no'
                        ## Get edit types ##
                        ##### This part is not relevant, vcf file has correct ref sequences already.... so c-U and A-I calls are what they look like without worrying about the ref nuc being from the wrong strand ####
                        if strand == '+':
                            if len(re) == len(alt1):
                                type1 = a1_pos_snp(re,alt1,type1)

                            if len(re) != len(alt1):
                                type1 = a1_pos_indel(re,alt1,type1)

                            if alt2 != None:
                                if len(re) == len(alt2):
                                    type2 = a2_pos_snp(re,alt2,type2)

                                if len(re) != len(alt2):
                                    type2 = a2_pos_indel(re,alt2,type2)
                            

                        if strand == '-':
                            if len(re) == len(alt1):
                                type1 = a1_neg_snp(re,alt1,type1)

                            if len(re) != len(alt1):
                                type1 = a1_neg_indel(re,alt1,type1)
                                
                            if alt2 != None:
                                if len(re) == len(alt2):
                                    type2 = a2_neg_snp(re,alt2,type2)

                                if len(re) != len(alt2):
                                    type2 = a2_neg_indel(re,alt2,type2)

                        rf.at[index,'string_id'] = string_id
                        rf.at[index,'ref_id'] = ref_id
                        rf.at[index,'tid'] = tid
                        rf.at[index,'exon_pos'] = exon_pos
                        rf.at[index,'trans_pos'] = trans_pos
                        rf.at[index,'single_iso'] = single_iso
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

