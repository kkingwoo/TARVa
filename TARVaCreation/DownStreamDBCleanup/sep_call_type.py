import os
import sqlite3 as sql


class CallTypes:
    @staticmethod
    def by_call_type(pos,dbp,gen,ad_pos_dict,cont_pos_dict,tidss):
        db_con = sql.connect(dbp,check_same_thread=False)
        db = db_con.cursor()
        ad_cts, cont_cts = 0,0
        ti_tupe = ''
        ftab = "fasta_tab"
        tids = f"SELECT * FROM {ftab} WHERE ref_id = ? AND tid = ? AND pos = ?"
        ti_list = []
        ti = db.execute(tids, (gen,tidss,pos)).fetchall()
        if ti:
            tis = ti[0][2]
            tp = ti[0][6]
            ti_tupe = (tis,tp)
        else:
            ti = 'na'
            tp = 'na'
            ti_tupe = (ti,tp)
        if gen in ad_pos_dict.keys():
            if pos in ad_pos_dict[gen].keys():
                ad_cts = ad_pos_dict[gen][pos]
        if gen in cont_pos_dict.keys():
            if pos in cont_pos_dict[gen].keys():
                cont_cts = cont_pos_dict[gen][pos]            
        db.close()
        db_con.close()
        return ti_tupe,ad_cts,cont_cts

    def check_ones_df(pos,gen,on,tidss,call,tps):
        unkown_list = []
        hgvs_list = []
        hgvs = 'na'
        out_list = []
        types = on[(on['position'] == pos) & (on['gene'] == gen)]
        alts = types['alt_seq'].tolist()
        refs = types['ref_seq'].tolist()
        uts = types['type'].tolist()
        rids = types['rid'].tolist()
        for alt,ref,ut,rid in zip(alts,refs,uts,rids):
            out_list1 = [gen,tidss,pos,alt,ref,ut,rid]
            out_row = [gen,tidss,pos,ut,call]
            out_list.append(out_row)
        tids = out_list1[1]
        if not tids == 'na':
            t = tidss.replace('_','.')
            ut = out_list1[5]
            ref,alt = out_list1[4],out_list1[3]
            rid = out_list1[6]
            if ut == 'other':
                if len(ref) == 1 and len(alt) == 1:
                    hgvs = f"{t}:c.{tps}{ref}>{alt}"
                    hgvs_list.append(hgvs)
                else:
                    unkown_out = [rid,gen,tidss,pos]
                    unkown_list.append(unkown_out)
            elif ut == 'insert':
                start,end = int(tps),int(tps)+len(alt)
                hgvs=f"{t}:c.{start}_{end}ins{alt}"
                hgvs_list.append(hgvs)
            elif ut == 'del':
                start,end = int(tps),int(tps)+len(ref)
                hgvs = f"{t}:c.{start}_{end}del"
                hgvs_list.append(hgvs)
            else:
                hgvs=f"{t}:c.{tps}{ref}>{alt}"
                hgvs_list.append(hgvs)
        unique_outs = list(set(map(tuple,out_list))) 
        out_list = [list(item) for item in unique_outs]

        return unkown_list,hgvs_list,out_list

    def novel_edits(pos,stra,chroms,chros,outs):
        re,da = '',''
        reds = red[(red['Position'].astype(int) == int(pos)) & (red['Strand'] == stra) & (red['Region'] == chroms)] #darns = darn[(darn['coordinate'].astype(int) == int(pos)) & (darn['strand'] == stra) & (darn['chrom'].astype(str) == str(chros))]
        if reds.empty:
            re = 'no'
        else:
             re = 'yes'

        if darns.empty:
            da = 'no'
        else:
            da = 'yes'
        
        outs.append(re)
        outs.append(da)

        return outs
