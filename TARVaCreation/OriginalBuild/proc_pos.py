import pandas as pd
import sqlite3 as sqlite

class ProcPos:
    @staticmethod
    def process_positions(gcounts,pcounts,dbp):
        pos_df = pd.DataFrame(columns=["ref_id","strand","pos","condition","ct","known"])
        db_con = sql.connect(dbp, check_same_thread=False)
        db = db_con.cursor()
        tab = "known_tab"
        chrom = ''
        for gs in gcounts:
            string_id,condition,rid,pos,a1_ad,a2_ad,dp,a1_edit_type,a2_edit_type,strand,chrom = gs[0],gs[1],gs[2],gs[3],gs[4],gs[5],gs[6],gs[7],gs[8],gs[9],gs[10]
            if not pos in pcounts[condition][strand]:
                pcounts[condition][strand][pos] = 0
            pcounts[condition][strand][pos] +=1
        
        checking_list = []
        for con in pcounts.keys():
            for stranding in pcounts[con].keys():
                for position in pcounts[con][stranding].keys(): 
                    pos_count = pcounts[con][stranding][position]
                    ck = f"SELECT COUNT(*) FROM {tab} WHERE strand = ? AND position = ?"
                    checking = db.execute(ck,(stranding,position)).fetchone()
                    checking_list.append(checking)

                            

        return checking_list                 #pos_df

