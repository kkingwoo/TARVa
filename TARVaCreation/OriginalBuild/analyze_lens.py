import sqlite3 as sql
from datetime import datetime
import numpy as np
from scipy.stats import mannwhitneyu

class GetLens:
    @staticmethod
    def from_gtf(db,ref,tid_list,gene_cts_list,trans_cts_list):
        state1 = f"Starting analysis for ref_id {ref} >>>> "+str(datetime.now())
        db_con = sql.connect(db,check_same_thread=False)
        dbs = db_con.cursor()
        gene_rids = {}
        gene_rids['Control'],gene_rids['AD'] = {},{}
        len_dict = {}
        gene_dict = {}
        gene_dict[ref] = {}
        gene_dict[ref]
        gene_dict[ref]['Positions'] = {}
        full_pos_list = []
        s_tab = 'sample_tab'
        g_tab = 'gtf_tab'
        samps_deets = []
        for i in range(0,len(gene_cts_list)):
            ad_samps_called,control_samps_called= 0,0
            tupe = gene_cts_list[i]
            posi,con,ct = tupe[0],tupe[1],tupe[2]
            from_samps = f"SELECT DISTINCT rid,a1_ad,a2_ad,dp FROM {s_tab} WHERE pos = ? AND condition = ? AND ref_id = ? GROUP BY rid"
            samps = dbs.execute(from_samps,(posi,con,ref)).fetchall()
            for s in samps:
                rid,a1_ad,a2_ad,dp = s[0],s[1],s[2],s[3]
                if not dp == None:
                    all_ad = a1_ad+a2_ad
                    if not posi in gene_dict[ref]['Positions'].keys():
                        gene_dict[ref]['Positions'][posi] = {}
                    if not con in gene_dict[ref]['Positions'][posi].keys():
                        gene_dict[ref]['Positions'][posi][con] = {}
                    if not rid in gene_dict[ref]['Positions'][posi][con].keys():
                        gene_dict[ref]['Positions'][posi][con][rid] = []
                        gene_rids[con][rid] = {}
                    out_tupe = (a1_ad,a2_ad,all_ad,dp)
                    gene_dict[ref]['Positions'][posi][con][rid].append(out_tupe)
                     
        state2 = f"Finished gene_dict for ref_id {ref} >>>> "+str(datetime.now())
        state3_list = []
        state4_list = []
        for t in range(0,len(trans_cts_list)):
            tid_rids = {}
            tupe = trans_cts_list[t]
            tid,posi,con,ct,wt = tupe[0],tupe[1],tupe[2],tupe[3],tupe[4]
            state3 = f"Starting trans_dict for ref_id {ref} : tid {tid} >>>>>>>"+str(datetime.now())
            state3_list.append(state3)
            lens_query = f"SELECT DISTINCT tid,trans_start,trans_end FROM {g_tab} WHERE tid= ? GROUP BY tid"
            lens = dbs.execute(lens_query, (tid,)).fetchall()
            for l in lens:
                start,end = l[1],l[2]
                full_pos_list.append(start)
                full_pos_list.append(end)
            state4 = f"Finished trans_dict for ref_id {ref} : tid {tid} >>>>>>>"+str(datetime.now())
            state4_list.append(state4)
            ##Do some stuff##
        
        
        full_pos_list.sort()
        full_set = list(set(full_pos_list))
        full_set.sort()

        g_start,g_end = full_set[0],full_set[-1]
        gene_len = g_end-g_start
        gene_pos_wt = 1/gene_len
        test_list = []
        for position in gene_dict[ref]['Positions'].keys():
            if int(position) in list(range(g_start,g_end)):
                for condition in gene_dict[ref]['Positions'][position].keys():
                    for rids in gene_dict[ref]['Positions'][position][condition].keys():
                        tupe_list = gene_dict[ref]['Positions'][position][condition][rids]
                        for tupe in tupe_list:
                            a1,a2,alls,dp = tupe[0],tupe[1],tupe[2],tupe[3]
                            if not 'a1' in gene_rids[condition][rids].keys():
                                gene_rids[condition][rids]['a1'],gene_rids[condition][rids]['a2'],gene_rids[condition][rids]['all'] = float(0.0),float(0.0),float(0.0)
                            a1s = float(a1/dp)*gene_pos_wt
                            a2s = float(a2/dp)*gene_pos_wt
                            allss = float(alls/dp)*gene_pos_wt
                            gene_rids[condition][rids]['a1']+=a1s
                            gene_rids[condition][rids]['a2']+=a2s
                            gene_rids[condition][rids]['all']+=allss

        ad,con = 29,30
        ad_all_list,con_all_list = [],[]
        ### This section is looking at total proportions per individual, per gene##
        #ad_a1_list,con_a1_list = [],[]
        #ad_a2_list,con_a2_list = [],[]
        cts_out = []
        for cons in gene_rids.keys():
            samp_ct=0
            for samp in gene_rids[cons].keys():
                samp_ct+=1
                all_tots = gene_rids[cons][samp]['all'] 
                #a1_tots = 
                #a2_tots = 
                if cons == "Control":
                    con_all_list.append(all_tots)
                if cons == "AD":
                    ad_all_list.append(all_tots)
            if cons == "Control":
                diff = con-samp_ct
                con_all_list.extend([float(0.0)] * diff)
                
            if cons == "AD":
                diff = ad-samp_ct
                ad_all_list.extend([float(0.0)] * diff)
                
        all_u,all_val = mannwhitneyu(ad_all_list, con_all_list, alternative='two-sided')
        
        state5 = f"Finished analysis for ref_id {ref} and all associated transcripts >>>> "+str(datetime.now())
        state6 = (ref,all_u,all_val)
        return state1,state2,state3_list,state4_list,state5,state6
