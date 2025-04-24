import sqlite3 as sql
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import pandas as pd
from datetime import datetime

class EachGene:
    @staticmethod
    def analyze_gene_info(dbp,gene):
        db_con = sql.connect(dbp,check_same_thread=False)
        db = db_con.cursor()
        
        pos_df = pd.DataFrame(columns=["ref_id","strand","pos","condition","ct","known"])
        s_tab = "sample_tab"
        conditions = ("Control","Control","AD","AD")
        types = ("A-I","C-U","A-I","C-U")
        ty = list(set(types))
        cons = ("Control","AD")
        strands = ('-','+')
        params = ','.join(['?'] * len(ty))
        gene_counts = f"SELECT rid,pos,strand,condition, COUNT(*) FROM {s_tab} WHERE string_id = ? AND ((a1_edit_type IN ({params})) OR (a2_edit_type IN ({params}))) GROUP BY condition,strand,pos"
        combined_vars = (gene,) + tuple(ty) + tuple(ty)
        gc = db.execute(gene_counts, combined_vars).fetchall()
        known = ''
        known_tab = "known_tab"
        for gs in gc:
            posing,stranding,conding,counting = gs[1],gs[2],gs[3],gs[4]
            ck = f"SELECT * FROM {known_tab} WHERE strand = ? AND position = ?"
            checking = db.execute(ck,(stranding,posing)).fetchone()
            if checking is not None:
                known = 'yes'
            else:
                known = 'no'
            out_list = [gene,stranding,posing,conding,counting,known]
            pos_df.loc[len(pos_df)] = out_list

        return pos_df
    

