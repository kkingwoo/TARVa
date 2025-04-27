import os
import sys
import csv
import sqlite3 as sql
import pandas as pd
import numpy as np
from dictionaries import Dictionaries
from datetime import datetime
from concurrent import futures
from make_sample_tabs import SampleTabs
from Bio import SeqIO
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
from raw_counts import RawCounts
from statsmodels.stats import multitest
from analyze_lens import GetLens

def create_fasta_table(fasta_file,dbp):
    print("Creating fasta_tab --> ",datetime.now())
    db_con = sql.connect(dbp,check_same_thread = False)
    db = db_con.cursor()
    gtf_table = 'gtf_tab'
    fasta_table = 'fasta_tab'
    tid_list = []
    tids = db.execute(f"SELECT DISTINCT tid FROM {gtf_table}").fetchall()
    for t in tids:
        ti = t[0]
        tid_list.append(ti)
    trans = set(tid_list)
    del tid_list
    db.execute("""
        CREATE TABLE IF NOT EXISTS %s (
            chrom varchar,
            pos integer,
            tid varchar,
            nuc varchar,
            strand varchar,
            exon_pos integer,
            trans_pos integer,
            exon_len integer,
            trans_len integer,
            per_exon_pos_wt float,
            per_trans_pos_wt float,
            exon varchar,
            string_id varchar,
            ref_id varchar,
            tid_pos varchar primary key);
        """ %(fasta_table))
    db_con.commit()
    
    key_list = ['chrom','pos','tid','nuc','strand','exon_pos','trans_pos','exon_len','trans_len','per_exon_pos_wt','per_trans_pos_wt','exon','string_id','ref_id','tid_pos']
    infast = open(fasta_file,'r')
    for record in SeqIO.parse(infast,'fasta'):
        insert_list = []
        all_exo_pos_list = []
        pos_list = []
        trans_pos=0
        ##seq_id is tid
        seq_id = record.id.split('|')[0].replace('.','_')
        if seq_id in trans:
            state = f"SELECT * FROM {gtf_table} WHERE tid = ?"
            get_info = db.execute(state,(seq_id,)).fetchall()
            sequence = str(record.seq)
            chrom,strand,exon,string_id,ref_id,nuc = '','','','','',''

            exon_start,exon_end,trans_start,trans_end = int(),int(),int(),int()

            for get in get_info:
                exon_position_list = []
                exon_count = 0
                chrom,ref_id,tid,strand,exon_start,exon_end,string_id,exon,trans_start,trans_end = get[0],get[1],get[2],get[3],get[4],get[5],get[6],get[7],get[8],get[9]
                for positions in range(exon_start,exon_end):
                    pos_list.append(positions)
                    exon_count+=1
                    exon_position_list.append(exon_count)
                exon_positions = tuple(exon_position_list)
                all_exo_pos_list.append(exon_positions)
            if len(pos_list) > 0:
                trans_len = len(pos_list)
                per_trans_pos_wt = float(1/trans_len)
                for e in range(0,len(all_exo_pos_list)):
                    exon = 'exon'+str(e+1)
                    exon_indices = all_exo_pos_list[e]
                    exon_len = len(exon_indices)
                    per_exon_pos_wt = float(1/exon_len)
                    for ex in range(0,exon_len):
                        if trans_pos == 0:
                            nuc = sequence[trans_pos]
                            trans_pos+=1
                        else:
                            if not trans_pos == 0:
                                try:
                                    nuc = sequence[trans_pos]
                                except IndexError:
                                    continue
                                trans_pos+=1
                        exon_pos = exon_indices[ex]
                        tp = trans_pos-1
                        pos = pos_list[tp]
                        tid_pos = tid+'_'+str(pos)
                        out_dict = {}

                        val_list = [chrom,pos,tid,nuc,strand,exon_pos,trans_pos,exon_len,trans_len,per_exon_pos_wt,per_trans_pos_wt,exon,string_id,ref_id,tid_pos]
                        insert_into = Dictionaries._loop_list_kvs(out_dict,key_list,val_list)

                        insert_list.append(insert_into)
            inserts = tuple(insert_list)
            fill_in = "INSERT OR IGNORE INTO %s VALUES(:chrom,:pos,:tid,:nuc,:strand,:exon_pos,:trans_pos,:exon_len,:trans_len,:per_exon_pos_wt,:per_trans_pos_wt,:exon,:string_id,:ref_id,:tid_pos)" %(fasta_table)
            db.executemany(fill_in,inserts)
            db_con.commit()
    db.execute(f"CREATE INDEX samp_parse ON {fasta_table}(chrom,pos,nuc);")
    db_con.commit()    
    print("FASTA table complete with index `samp_parse` on chrom,pos,nuc --> ",datetime.now())
    return 

def make_gtf_table(gtf_file,dbp):
    print("Starting gtf_tab -->",datetime.now())
    insert_list = []
    iso_dict = {}
    ## 1.) Create db connection and assign table name
    db_con = sql.connect(dbp)
    db = db_con.cursor()
    ## 2.) Create table
    table_name = 'gtf_tab'
    #db.execute(f"DROP TABLE IF EXISTS {table_name}")
    db.execute("""
        CREATE TABLE IF NOT EXISTS %s (
            chrom varchar,
            ref_id varchar,
            tid varchar,
            strand varchar,
            exon_start integer,
            exon_end integer,
            string_id varchar,
            exon varchar,
            trans_start integer,
            trans_end integer,
            tid_exon varchar primary key);
        """ %(table_name))
    db_con.commit()
    key_list = ['chrom','ref_id','tid','strand','exon_start','exon_end','string_id','exon','trans_start','trans_end','tid_exon']
    insert_list = []
    gtf = open(gtf_file, 'r')
    ref_id_dict = {}
    for l in gtf.readlines(): 
        if l.startswith('#'):
            continue
        else:
            l = l.split()
            chro = l[0]
            if not chro in iso_dict.keys():
                iso_dict[chro] = {}
            feature = l[2]
            tid = l[11].split('"')[1].replace('.','_')
            string_id = l[9].split('"')[1].replace('.','_')
            if not string_id in ref_id_dict.keys():
                ref_id_dict[string_id] = {}
                ref_id_dict[string_id][tid] = {}
            if not tid in ref_id_dict[string_id].keys():
                ref_id_dict[string_id][tid] = {}
            if feature == "transcript":
                trans_start = int(l[3])
                trans_end = int(l[4])+1
                trans_tupe = (trans_start,trans_end)
                ref_id = l[-1].split('"')[1].replace('.','_')
                ref_id_dict[string_id][tid][ref_id] = trans_tupe
                if not string_id in iso_dict[chro].keys():
                    iso_dict[chro][string_id] = 1
                else:
                    iso_dict[chro][string_id]+=1
            if feature == "exon":
                exon_num = 'exon'+l[13].split('"')[1]
                ref_id = list(ref_id_dict[string_id][tid].keys())[0]
                out_dict = {}
                exon_start = int(l[3])
                exon_end = int(l[4])+1
                strand = l[6]
                tid_exon = tid+'_'+exon_num
                trans_start,trans_end = ref_id_dict[string_id][tid][ref_id][0],ref_id_dict[string_id][tid][ref_id][1]
                val_list = [chro,ref_id,tid,strand,exon_start,exon_end,string_id,exon_num,trans_start,trans_end,tid_exon]
                insert_into = Dictionaries._loop_list_kvs(out_dict,key_list,val_list)
                insert_list.append(insert_into)

    inserts = tuple(insert_list)
    fill_in = "INSERT OR IGNORE INTO %s VALUES(:chrom,:ref_id,:tid,:strand,:exon_start,:exon_end,:string_id,:exon,:trans_start,:trans_end,:tid_exon)" %(table_name)
    db.executemany(fill_in,inserts)
    db_con.commit()

    print("GTF table complete -->",datetime.now())
    return iso_dict


def get_matched_files(wp,rp):
    match_dict  = {}
    match_dict['AD'],match_dict['Control'] = {},{}
    rna_ad,rna_con = rp+'AD/',rp+'Control/'
    rna_list = [rna_ad,rna_con]
    wgs_ad,wgs_con = wp+'AD/',wp+'Control/'
    wgs_list = [wgs_ad,wgs_con]
    
    for r in range(0,len(rna_list)):
        for files in os.listdir(rna_list[r]):
            rid = files.split('.')[0]
            zf = rna_list[r]+files
            if r == 0:
                match_dict["AD"][rid] = []
                match_dict["AD"][rid].append(zf)
            if r == 1:
                match_dict["Control"][rid] = []
                match_dict["Control"][rid].append(zf)
    
    for w in range(0,len(wgs_list)):
        for vcfs in os.listdir(wgs_list[w]):
            rids = vcfs.split('_')[0] 
            if w == 0:
                if rids in match_dict['AD'].keys():
                    subs = wgs_list[w]+vcfs+'/'
                    for files in os.listdir(subs):
                        if files.endswith("contaminates_filtered.vcf"):
                            vf  = subs+files
                            match_dict["AD"][rids].append(vf)                 
        
            if w == 1:
                if rids in match_dict["Control"].keys():
                    subs = wgs_list[w]+vcfs+'/'
                    for files in os.listdir(subs):
                        if files.endswith("contaminates_filtered.vcf"):
                            vf  = subs+files
                            match_dict["Control"][rids].append(vf)
    return match_dict

def make_samp_tabs(mf,gdict,db_path):

    samp_tab = 'sample_tab'
    for c in mf.keys():
        final_results = []
        with futures.ProcessPoolExecutor(max_workers=32) as mst:
            print('************Creating sample table for ',c,' group --> ',datetime.now())
            wait_for = [mst.submit(SampleTabs.sample_tabs,c,s,mf,db_path,gdict) for s in mf[c].keys() ]
            for fu in futures.as_completed(wait_for):
                current = fu.result()
                final_results.append(current)            
        db_con = sql.connect(db_path,check_same_thread=False)
        for finals in final_results:
            finals.to_sql(samp_tab,db_con,index=False,if_exists='append')
    print('SAMPLE table has been created and populated\n')
                     
    return


def make_gene_list(db_path):
    
    gtfs = 'gtf_tab'
    db_con = sql.connect(db_path,check_same_thread=False)
    db = db_con.cursor()
    genes = db.execute(f"SELECT DISTINCT ref_id FROM {gtfs}").fetchall()
    gene_list = [g[0] for g in genes]
    db.close()
    db_con.close()
    
    return gene_list

def make_known_dfs(db_path):

    known_tab = "known_tab"
    db_con  = sql.connect(db_path,check_same_thread=False)
    known_edits = open('TABLE1_hg38.txt','r')
    known = csv.reader(known_edits,delimiter='\t')
    next(known)
    known_list = []
    for k in known:
        position = k[1]
        strand = k[4]
        known_list.append([str(strand),str(position)])
    known_df = pd.DataFrame(known_list, columns=['strand','position'])
    known_df.to_sql(known_tab,db_con,index=False,if_exists='append')
    db_con.close()
    return                       

def get_rids_tids(db_path,g_list):

    samp_tab = "sample_tab"
    db_con = sql.connect(db_path,check_same_thread=False)
    db = db_con.cursor()
    gene_params = ','.join(['?'] * len(g_list))
    tids = f"SELECT DISTINCT ref_id,tid FROM {samp_tab} WHERE ref_id IN ({gene_params})"
    ti = db.execute(tids, g_list).fetchall()
    del g_list
    del gene_params
    rids = f"SELECT DISTINCT rid FROM {samp_tab} WHERE condition = ?"
    get_ads = db.execute(rids,("AD",)).fetchall()
    a_rids = [a[0] for a in get_ads]
    get_controls = db.execute(rids,("Control",)).fetchall()
    c_rids = [co[0] for co in get_controls]
    id_tupe = (c_rids,a_rids)
    gene_dict = {}
    for tu in ti:
        key = tu[0]
        val = tu[1]
        if not key in gene_dict.keys():
            gene_dict[key] = []
        gene_dict[key].append(val)
        
    db.close()
    db_con.close()

    return id_tupe,gene_dict

def make_per_pos_gene_file(db_path,ids,g_dict):
    positions_header = ["strand","position","condition","ct","known"]
    with open('per_position_gene.csv','w') as position_file:
        pos_writer = csv.writer(position_file)
        pos_writer.writerow(positions_header)
        position_file.close()
    pos_append = open('per_position_gene.csv','a')
    out_write = csv.writer(pos_append)
    samp_tab = "sample_tab"
    known_tab = "known_tab"
    gtfs = 'gtf_tab'
    db_con = sql.connect(db_path,check_same_thread=False)
    db = db_con.cursor()
    conditions = ("Control","Control","AD","AD")
    types = ("A-I","C-U","A-I","C-U")
    ty = list(set(types))
    cons = ("Control","AD")
    strands = ('-','+')
    concat_list = []
    params = ','.join(['?'] * len(ty))
    genes = db.execute(f"SELECT condition,pos,strand, COUNT(*) AS count FROM {samp_tab} GROUP BY condition,pos,strand").fetchall()
    for ge in genes:
        known = ''
        cond,po,stra,ct = ge[0],ge[1],ge[2],ge[3]
        ck = f"SELECT * FROM {known_tab} WHERE strand = ? AND position = ?"
        checking = db.execute(ck,(stra,po)).fetchone()
        if checking is not None:
            known = 'yes'
        else:
            known = 'no'
        out_list = [stra,po,cond,ct,known]
        out_write.writerow(out_list)

    db.close()
    db_con.close()

    return

def get_fast_inf(dbp,cf):
    iso_dict = {}
    trans_pos_cts = {}
    gene_pos_cts = {}
    in_file = open(cf,'r')
    in_reader = csv.reader(in_file)
    next(in_reader)
    db_con = sql.connect(dbp,check_same_thread=False)
    db = db_con.cursor()
    tab_name = 'fasta_tab'
    print("starting fasta_tab row parsing >>>> ",datetime.now())
    for i in in_reader:
        stra,posi,con,ct = i[0],i[1],i[2],i[3]
        query = f"""
        SELECT DISTINCT ref_id,tid,per_trans_pos_wt
        FROM {tab_name}
        WHERE strand= ? AND pos = ? GROUP BY ref_id"""
        db_res = db.execute(query, (stra, posi)).fetchall()
        for r in db_res:
            refs,tids,wts = r[0],r[1],r[2]
            gene_tupe,trans_tupe = (posi,con,ct,wts),(tids,posi,con,ct,wts)
            if not refs in iso_dict.keys():
                iso_dict[refs] = []
                gene_pos_cts[refs] = []
                trans_pos_cts[refs] = []
            if not gene_tupe in gene_pos_cts[refs]:
                gene_pos_cts[refs].append(gene_tupe)
            if not trans_tupe in trans_pos_cts[refs]:
                trans_pos_cts[refs].append(trans_tupe)
            if not tids in iso_dict[refs]:
                iso_dict[refs].append(tids)
    
    print("fasta_tab row parsing finished >>>>> ",datetime.now())
    db.close()
    db_con.close()
    return iso_dict,gene_pos_cts,trans_pos_cts

def per_pos_analyze(in_dat,dbp,tops):
    isos,gene_pos,trans_pos = in_dat[0],in_dat[1],in_dat[2]
    db_con = sql.connect(dbp)
    dbs = db_con.cursor()
    len_dict = {}
    gene_lens = {}
    outs = []
    process_num = 0
    with futures.ProcessPoolExecutor(max_workers=32) as mst:
        print('************Getting lengths for per_position analysis --> ',datetime.now())
        wait_for = [mst.submit(GetLens.from_gtf,dbp,ref_id,isos[ref_id],gene_pos[ref_id],trans_pos[ref_id]) for ref_id in tops if ref_id in isos.keys() and ref_id in gene_pos.keys() and ref_id in trans_pos.keys()]
        for fu in futures.as_completed(wait_for):
            current = fu.result()
            print(current[0],'\n','\n')
            print(current[1],'\n','\n')
            for s in range(0,len(current[2])):
                print(current[2][s],'\n')
                print(current[3][s],'\n','\n')
            print(current[4])
            outs.append(current[5])
            
    print("Per position analysis is now complete >>>>>", datetime.now())
    
    with open('all_proportions_pvals.csv','w') as props_vals:
        props_writer = csv.writer(props_vals)
        head = ['gene_id','u_stat','p_val']
        props_writer.writerow(head)
        for o in outs:
            ref_id,u_stat,p_val = o[0],o[1],o[2]
            outrow = [ref_id,u_stat,p_val]
            props_writer.writerow(outrow)
    return
    

if __name__=='__main__':
    gtf_file_path = sys.argv[1] 
    fasta_path = sys.argv[2] 
    db_path = sys.argv[3]
    genes = sys.argv[4]
    isos = sys.argv[5]
    wgs_path = sys.argv[6]
    rna_path = sys.argv[7]
    string_path = sys.argv[8]
    counts_file = sys.argv[9]
    gtfs = make_gtf_table(gtf_file_path,db_path)
    fastas = create_fasta_table(fasta_path,db_path)
    g_tables = create_gene_tables(gtfs,db_path,genes_of_interest)
    length_dictionary = len_dict(gtfs)
    matched_files = get_matched_files(wgs_path,rna_path)
    samp_tabs = make_samp_tabs(matched_files,gtfs,db_path)
    genes = make_gene_list(db_path)
    all_ids = get_rids_tids(db_path,genes)
    k_dfs = make_known_dfs(db_path)
    pos_file = make_per_pos_gene_file(db_path,all_ids[0],all_ids[1])
    fast_inf = get_fast_inf(db_path,counts_file)
    gtf_inf = per_pos_analyze(fast_inf,db_path,genes)
    
