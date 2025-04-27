import csv
import os
import sys
from datetime import datetime
from concurrent import futures
from type_by_gene import EditTypes
import scipy.stats as stats
import sqlite3 as sql
import pandas as pd
from adjust_lists import ListLens
from scipy.stats import mannwhitneyu
from sep_call_type import CallTypes
from parse_vep import ParseVep
import pandas as pd

def make_new_UUIDs(d):
    id_dict = {}
    id_dict['AD'] = {}
    id_dict['Control'] = {}
    db_con = sql.connect(d)
    db = db_con.cursor()
    stab = "sample_tab"
    old_ids = db.execute(f"SELECT DISTINCT rid,condition FROM {stab} GROUP BY rid").fetchall()
    act,cct = 1,1
    for old in old_ids:
        ori = old[0]
        condition = old[1]
        if condition =="AD":
            new_id= condition+'_'+str(act)
            id_dict['AD'][ori] = new_id
            act+=1
        if condition == "Control":
            new_id= condition+'_'+str(cct)
            id_dict["Control"][ori] = new_id
            cct+=1

    return id_dict

def per_genes(d,i,o1,o2,newids,red):
    out = []
    genes = []
    reds = pd.read_csv(red,low_memory=False,delimiter='\t')
    
    orig_clin = pd.read_csv("/scratch/kkingwoo/ADProj/ROSMAP_BulkBrain_2PassMapped/ROSMAP_All_Metadata/ROSMAP_clinical.csv")
    ins = open(i,'r')
    in_read = csv.reader(ins)
    next(in_read)
    for re in in_read:
        pval = float(re[3])
        if pval <= float(0.05):
            gene = re[0].replace('.','_') 
            genes.append(gene)
    with futures.ProcessPoolExecutor(max_workers=32) as mst:
        print('******************Starting edit-type analysis for genes with significant differences in editing levels between the two conditions -->',datetime.now())
        wait_for = [mst.submit(EditTypes.edits,d,g) for g in genes]
        for fu in futures.as_completed(wait_for):
            current = fu.result()
            out.append(current)
    
    for ou in out:
        print(ou)
    
    ##unique_outs,all_info_outs = [],[]
    #ensg_list,mstrg_list = [],[]
        
    hugo_reading = pd.read_csv('top_gene_names_for_Enrichr.txt')

    #all_info_df = pd.read_csv('Main_All_Info.csv')
    #all_info_df = all_info_df.drop('p-value',axis=1)
    

    dist_type_list = ['AD_only','Control_only','Primarily_AD','Primarily_Control','Both']
    type_file_list = []
    for di in range(0,len(dist_type_list)):
        pref = dist_type_list[di]
        fil = f"{pref}_For_Enrichr.txt"
        type_file_list.append(fil)
    ad_only_list,control_only_list,prim_ad_list,prim_con_list,both_list = [],[],[],[],[]
    type_lists_list = [ad_only_list,control_only_list,prim_ad_list,prim_con_list,both_list]


    #hugo_out = open('refids_for_HUGO.csv','w')
    #hugo_write = csv.writer(hugo_out)

    #hgvs_out = open("All_HGVS_for_VEP.txt",'w')
    #hgvs_write=  csv.writer(hgvs_out,delimiter='\t')
    
    distinct_outs = []
    out_string_list = []
    #write_header = True
    #cols = ['con','rids','poss','refid','ref_base','a1ad','dp','ref','alt1','a1type','strand','chrom','hgvs','a1prop','edit_dist']
    #all_info_df = pd.DataFrame(columns=cols)
    conn = sql.connect(d)
    table = 'All_Info_TopGenes_tab'
    co = conn.cursor()
    co.execute(f"DROP TABLE IF EXISTS {table}")
    conn.commit()
    for outing in out:
        refid = outing[0]
        gene_dist_type = outing[1]
        out_stats = outing[2]
        unique_positions_dict = outing[3]
        check = outing[4]
        hgvs_list = outing[5]
        out_string = outing[6]
        out_string_list.append(out_string)
        dist_dict = outing[7]
        #hugo_outs.append(refid)
     
        unique_outs = []
        #write_header = True

        cols = ['con','rids','poss','refid','ref_base','a1ad','dp','ref','alt1','a1type','strand','chrom','hgvs','a1prop','edit_dist']
        all_info_df = pd.DataFrame(columns=cols)

        #all_info_df = pd.read_csv('Main_All_Info.csv')
        for ch in check:
            ch = list(ch)
            ch.append(gene_dist_type)
            all_info_df.loc[len(all_info_df)] = ch

        if out_stats:
            all_info_df['p_value'] = None
            for ou in out_stats:
                p_val,rd,t = ou[3],ou[0],ou[1]
                mask = (all_info_df['refid'] == rd) & (all_info_df['a1type'] == t)
                all_info_df.loc[mask, 'p_value'] = p_val

        #for hg in hgvs_list:
        #    all_hgvs_list.append(hg)
        
        all_info_df['redi']=None
        for strands in unique_positions_dict.keys():
            if unique_positions_dict[strands].keys():
                for chrom in unique_positions_dict[strands].keys():
                    for po in unique_positions_dict[strands][chrom]:
                        red = reds[(reds['Position'].astype(int) == int(po)) & (reds['Strand'] == strands) & (reds['Region'] == chrom)]
                        re = 'yes' if not red.empty else 'no'
                        mask = (all_info_df['strand'] == strands) & (all_info_df['chrom'] == chrom) & (all_info_df['poss'].astype(int) == int(po))
                        all_info_df.loc[mask,'redi'] = re
        print(all_info_df)
            
        rids_list = all_info_df['rids'].unique().tolist()
        all_info_df['msex'],all_info_df['apoe_genotype'],all_info_df['age_death'],all_info_df['braaksc'],all_info_df['dcfdx_lv'], all_info_df['ceradsc'] = None,None,None,None,None,None 
        for rid in rids_list:
            clin_data = orig_clin[orig_clin['individualID'] == rid][['msex','apoe_genotype','age_death','braaksc','dcfdx_lv','ceradsc']].values.tolist()
            if clin_data:
                mask = all_info_df['rids'] == rid
                all_info_df.loc[mask, ['msex', 'apoe_genotype', 'age_death', 'braaksc', 'ceradsc','dcfdx_lv']] = clin_data[0]

        all_info_df['gene_name'] = None

        ensg_list = all_info_df['ref_base'].unique().tolist()
        for ensg in ensg_list:
            hugs = hugo_reading[hugo_reading['Gene_stable_ID'] == ensg]['Gene_name'].values.tolist()
            if hugs:
                gene_name = hugs[0]
                mask = all_info_df['ref_base'] == ensg
                all_info_df.loc[mask,'gene_name'] = gene_name

            dist = all_info_df[all_info_df['ref_base']== ensg]['edit_dist'].values.tolist()
            d = dist[0] if dist else None
            if d in dist_type_list:
                ty = dist_type_list.index(d)
                the_list = type_lists_list[ty]
                the_list.append(gene_name)
            
        all_info_df.to_sql('All_Info_TopGenes_tab',conn,if_exists='append',index=False)
        print(f"All info for ensembl gene {refid} complete >>>>",datetime.now()) 
        
        ## Create a csv frile from the current dataframe, read in for the next testing step and just continue to add to it:

        #all_info_df.to_csv("Main_All.csv", mode='a',header=write_header, index=False)
        #write_header=False
        
    conn.commit()
    conn.close()
    #query = co.execute(f"SELECT DISTINCT edit_dist,gene_name FROM {table} GROUP BY edit_dist,gene_name").fetchall()
    #for qu in query:
    #    typ = qu[0]
    #    gn = qu[1]
    #    if not gn is None:
    #        loc = dist_type_list.index(typ)
    #        li = type_lists_list[loc]
    #        li.append(gn)
    #for index in range(0,len(type_lists_list)):
    #    fi = type_file_list[index]
    #    fil = open(fi,'w')
    #    fi_writer = csv.writer(fil,delimiter = '\t')
    #    lis = type_lists_list[index]
    #    for g in lis:
    #        fi_writer.writerow([g])
    

    #for types in range(0,len(dist_type_list)):
    #    out_list,out_file = type_lists_list[types],type_file_list[types]
    #    writing =  open(out_file,'w')
    #    for h in out_list:
    #      writing.write(h)
        
    #    if hugos.startswith("ENSG"):
    #        hugs = hugos.split('_')[0]
    #        hugo_write.writerow([hugs])
    #    if hugos.startswith("MS"):
    #        print(f"{hugos} is stringtie-identified as novel and thus has no associated HUGO id\n")
    #for ids in all_hgvs_list:
    #    hgvs_write.writerow([ids])
    
    return 

if __name__=='__main__':
    dbp = sys.argv[1]
    in_file = sys.argv[2]
    out_file = sys.argv[3]
    out_file2 = sys.argv[4]
    redi = sys.argv[5]
    uuid_dict = make_new_UUIDs(dbp)
    types_analyzed = per_genes(dbp,in_file,out_file,out_file2,uuid_dict,redi)
