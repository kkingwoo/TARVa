import sqlite3 as sql


class GetCounts:
    @staticmethod
    def get_counts(db,se,ct,tab):
        db_con = sql.connect(db,check_same_thread = False)
        d = db_con.cursor() 
        out_list = []
        ct+=1
        c = str(ct)
        ensg,poss = se[0],se[1]
        mod = f"mod_{c}"
        rids2 = f"""SELECT tissue,condition,ref_id,pos,COUNT(*) FROM {tab} WHERE ref_id = ? AND pos = ? GROUP BY tissue,condition,ref_id,pos"""
        r2 = d.execute(rids2,(ensg,poss)).fetchall()
        if r2:
            for twos in r2:
                out_list.append(mod)
                out_list.append(twos)
        
        return out_list
