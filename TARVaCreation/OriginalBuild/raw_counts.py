import sqlite3 as sql
import numpy as np
from scipy.stats import ttest_ind
from datetime import datetime
from statsmodels.stats import multitest
class RawCounts:
    @staticmethod
    def raw_counts_per_gene(dbp,ge,params,a,c,s_tab,ty):
        outrow = []
        if not type(ge) == None:
            db_con = sql.connect(dbp,check_same_thread=False)
            db = db_con.cursor()
            ad_ai,ad_cu,co_ai,co_cu = [],[],[],[]
            infos = f"SELECT condition,rid,a1_edit_type,a2_edit_type,COUNT(*) FROM {s_tab} WHERE ref_id = ? AND (a1_edit_type IN ({params}) OR a2_edit_type IN ({params})) GROUP BY ref_id,condition,rid,a1_edit_type,a2_edit_type"
            inf = db.execute(infos,(ge,)+tuple(ty*2)).fetchall()
            for i in inf:
                con,a1_edit_type,a2_edit_type,ct = i[0],i[2],i[3],i[4]
           ##Create module for this?
                if con == "AD":
                    if a1_edit_type == 'A-I' or a2_edit_type == 'A-I':
                        ad_ai.append(ct)
                    if a1_edit_type == 'C-U' or a2_edit_type == 'C-U':
                        ad_cu.append(ct)
                if con == "Control":
                    if a1_edit_type == 'A-I' or a2_edit_type == 'A-I':
                        co_ai.append(ct)
                    if a1_edit_type == 'C-U' or a2_edit_type == 'C-U':
                        co_cu.append(ct)
            #outrow = []
            len_adai,len_adcu,len_coai,len_cocu = str(len(ad_ai)),str(len(ad_cu)), str(len(co_ai)),str(len(co_cu))
            if  int(len_adai) >= 10 or int(len_coai) >= 10 or int(len_adcu) >= 10 or int(len_cocu) >= 10:

                adai_diff = a - int(len_adai)
                adcu_diff = a - int(len_adcu)
                coai_diff = c - int(len_coai)
                cocu_diff = c - int(len_cocu)
                for aaid in range(0,adai_diff):
                    ad_ai.append(0)
                for acud in range(0,adcu_diff):
                    ad_cu.append(0)
                for caid in range(0,coai_diff):
                    co_ai.append(0)
                for coud in range(0,cocu_diff):
                    co_cu.append(0)
        

                outrow = []
                ai_stat,ai_pval = ttest_ind(ad_ai,co_ai,equal_var=False)
                cu_stat,cu_pval = ttest_ind(ad_cu,co_cu,equal_var=False)
                out_tupe = (ge,ai_pval,len_adai,len_coai,cu_pval,len_adcu,len_cocu)
                outrow.append(out_tupe)
        
        return outrow
    

