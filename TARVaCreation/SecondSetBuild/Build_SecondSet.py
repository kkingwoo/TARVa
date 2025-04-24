import os
import sys
import csv
import json
import sqlite3 as sql
import pandas as pd
import numpy as np
import math
from scipy.stats import fisher_exact
from dictionaries import Dictionaries
from datetime import datetime
from concurrent import futures
from make_secondset_tabs import SampleTabs
from statsmodels.stats.multitest import multipletests
from get_second_set_counts import GetCounts

def make_modid_dict(firsts,ssmod):
    modid_dict,dump_dict = {},{}
    f = pd.read_csv(firsts)
    uniques = f[['ensg', 'pos']].drop_duplicates().to_records(index=False)
    unique_pairs_list = [tuple(row) for row in uniques]
    ct = 0
    for tu in unique_pairs_list:
        ct+=1
        num = str(ct)
        modid = f"mod{num}"
        ensg,po = tu[0],tu[1]
        modid_dict[modid] = {ensg:po}
        dump_dict[modid] = {ensg:str(po)}
        with open(ssmod,"w") as json_file:
            json.dump(dump_dict,json_file)
    return modid_dict

def make_samp_tabs(db,filepaths,tissues,clindats,out_file,md):
    patt = 'filtered.vcf'
    samp_tab = 'SecondSet_AllInfo_tab'
    final_results = []
    for i in range(0,len(filepaths)):
        vpa,tiss = filepaths[i],tissues[i]
        vfiles = [f for f in os.listdir(vpa) if f.endswith(patt)]
        with futures.ProcessPoolExecutor(max_workers=16) as mst:
            print('************Creating sample table for ',tiss,' group --> ',datetime.now())
            wait_for = [mst.submit(SampleTabs.sample_tabs,vpa,fi,tiss,clindats,db) for fi in vfiles]
            for fu in futures.as_completed(wait_for):
                current = fu.result()
                final_results.append(current)            
    db_con = sql.connect(db)
    dbcursor= db_con.cursor()
    #db_con = sql.connect(db_path,check_same_thread=False)
    for finals in final_results:
        finals.to_sql(samp_tab,db_con,index=False,if_exists='append')
    print('SAMPLE table for second set has been created and populated\n')
    
    return

def all_tissues_conditions(db,outfile,mdi):
    db_con =  sql.connect(db)
    db_cursor = db_con.cursor()
    tissue = 'bulkbrain'
    first_set,second_set = 'All_Info_TopGenes_tab','SecondSet_AllInfo_tab'
    
    key_vals_main = [
            (main_key, key, val)
            for main_key, sub_dict in mdi.items()
            for key, val in sub_dict.items()
         ]
    keys_vals = [(main_key,ref_base, poss) for main_key, ref_base, poss in key_vals_main]
    run_queries = []
    for tupe in keys_vals:
        mains,ensg,poss = tupe
        poss = int(poss)
        
        query_one = f"SELECT ref_base,poss,con,rids FROM {first_set} WHERE ref_base= ? AND poss = ?"
        run_query1 = db_cursor.execute(query_one,(ensg,poss)).fetchall()
        
        query_two = f"SELECT ref_id,pos,tissue,condition,rid FROM {second_set} WHERE ref_id = ? AND pos = ?"
        run_query2 = db_cursor.execute(query_two,(ensg,poss))
        if run_query1:
            for run1 in run_query1:
                outs = [run1[0],run1[1],tissue,run1[2],run1[3],mains]
                run_queries.append(outs)
        if run_query2:
            for run2 in run_query2:
                outs = [run2[0],run2[1],run2[2],run2[3],run2[4],mains]
                run_queries.append(outs)

    cols = ['ref','pos','tissue','condition','rid','main_key']
   
    df = pd.DataFrame(run_queries, columns = cols)
    df.to_csv(outfile,mode='a',header=False,index=False)

    return

def get_sample_sizes(db):
    dbcon = sql.connect(db)
    cursor = dbcon.cursor()
    sample_tab = 'SecondSet_AllInfo_tab'
    query = f"""SELECT condition,tissue, COUNT(DISTINCT rid) FROM {sample_tab} GROUP BY condition,tissue"""
    que = cursor.execute(query).fetchall()
    for q in que:
        print(q)

    return

def run_analysis_for_volcano(ofile):
    o = pd.read_csv(ofile)
    counts_dict = {}
    sample_counts = {"bulkbrain_AD":29,"bulkbrain_Control":30,"bulkbrain_MCI":27,"monocyte_AD":34,"monocyte_Control":89,"monocyte_MCI":38}
    one,two,three,four,five,six = 'bulkbrain_AD','bulkbrain_Control','bulkbrain_MCI','monocyte_AD','monocyte_Control','monocyte_MCI'
    comparison_tupes = [(one,four),(one,five),(one,six),(two,four),(two,five),(two,six),(three,four),(three,five),(three,six)]
    counts_dict["mods"] = []
    con_list,tissue_list = ['AD','Control','MCI'],['bulkbrain','monocyte']
    unique_vals = pd.unique(o['mod_id'].tolist())
    for u in unique_vals:
        counts_dict["mods"].append(u)
        for ti in tissue_list:
            for co in con_list:
                key_name = f"{ti}_{co}"
                if not key_name in counts_dict.keys():
                    counts_dict[key_name] = []
                flt = o[ 
                        (o['tissue'] == ti) &
                        (o['condition'] == co) & 
                        (o['mod_id'] == u)
                        ]
                
                counts_dict[key_name].append(len(flt))
    counting = pd.DataFrame(counts_dict)
    
    results = []
    for mod in counting["mods"]:
        mod_results = {"mod_id":mod}
        for cond1,cond2 in comparison_tupes:
            if cond1 in counting and cond2 in counting:
                count1  =  counting.loc[counting["mods"] == mod, cond1].values
                count2  =  counting.loc[counting["mods"] == mod, cond2].values

                if len(count1) == 0 or len(count2) == 0:
                    continue
                total1 = sample_counts[cond1]
                total2 = sample_counts[cond2]
                
                table = [[count1[0], total1 - count1[0]], [count2[0], total2 - count2[0]]]
                oddsratio, pvalue = fisher_exact(table)

                mod_results[f"{cond1}_vs_{cond2}_OddsRatio"] = oddsratio
                mod_results[f"{cond1}_vs_{cond2}_P_Value"] = pvalue
        results.append(mod_results)
    results_df = pd.DataFrame(results)
    
    if not results_df.empty:
        p_value_columns = [col for col in results_df.columns if "P_Value" in col]
        
        for col in p_value_columns:
            adjusted_pvals = multipletests(results_df[col], method="fdr_bh")[1]
            adjusted_pvals = [1e-300 if math.isinf(p) else p for p in adjusted_pvals]
            results_df[f"{col}_Adjusted"] = adjusted_pvals
    results_df.to_csv("fisher_results.csv",index=False)
    
    print(results_df)

    return 

def get_samp_info(db,clinical):
    info_dict = {}
    missing_vals_dict = {}
    clin = pd.read_csv(clinical)
    #samp_dat = clin_info.loc[clin_info['individualID'] == samp,["msex","age_death","apoe_genotype","braaksc","ceradsc","dcfdx_lv"]].iloc[0].tolist()
    bb_set,second_set = "All_Info_TopGenes_tab","SecondSet_AllInfo_tab"
    bb_con,ss_con = sql.connect(db,check_same_thread = False),sql.connect(db,check_same_thread = False)
    bb_cursor,ss_cursor = bb_con.cursor(),ss_con.cursor()
    bq,sq = f"SELECT DISTINCT rids,con,age_death,apoe_genotype,braaksc,ceradsc,msex FROM {bb_set} GROUP BY rids",f"SELECT DISTINCT rid,condition,age,apoe,braak,cerad,sex,tissue FROM {second_set} GROUP BY rid"
    b_sub,q_sub = bb_cursor.execute(bq), ss_cursor.execute(sq)
    sub_keys = ['sex','age','apoe','braak','cerad']
    for rids in b_sub:
        condition,age,apoe,braak,cerad,sex = rids[1],rids[2],rids[3],rids[4],rids[5],rids[6]
        tissue = "bulkbrain"
        if not condition in info_dict:
            info_dict[condition] = {}
            missing_vals_dict[condition] = {}
        if not tissue in info_dict[condition]:
            missing_vals_dict[condition][tissue] = {}
            info_dict[condition][tissue] = {}
            info_dict[condition][tissue]['90+'] = {}
            info_dict[condition][tissue]['under_90'] = {}
            for i in sub_keys:
                info_dict[condition][tissue]['90+'][i] = []
                info_dict[condition][tissue]['under_90'][i] = []
                missing_vals_dict[condition][tissue][i] = 0
        if age == '90+':
            subs = '90+'
            if sex == None:
                missing_vals_dict[condition][tissue]['sex']+=1
            else:
                info_dict[condition][tissue][subs]['sex'].append(sex)
            info_dict[condition][tissue][subs]['age'].append(age)
            if apoe == None:
                missing_vals_dict[condition][tissue]['apoe']+=1
            else:
                info_dict[condition][tissue][subs]['apoe'].append(apoe)
            if braak == None:
                missing_vals_dict[condition][tissue]['braak']+=1
            else:
                info_dict[condition][tissue][subs]['braak'].append(braak)
            if cerad == None:
                missing_vals_dict[condition][tissue]['cerad']+=1
            else:
                info_dict[condition][tissue][subs]['cerad'].append(cerad)
        if age == None:
            missing_vals_dict[condition][tissue]['age']+=1
        if age != '90+' and age != None:
            subs = 'under_90'
            if sex == None:
                missing_vals_dict[condition][tissue]['sex']+=1
            else:
                info_dict[condition][tissue][subs]['sex'].append(sex)
            info_dict[condition][tissue][subs]['age'].append(age)
            if apoe == None:
                missing_vals_dict[condition][tissue]['apoe']+=1
            else:
                info_dict[condition][tissue][subs]['apoe'].append(apoe)
            if braak == None:
                missing_vals_dict[condition][tissue]['braak']+=1
            else:
                info_dict[condition][tissue][subs]['braak'].append(braak)
            if cerad == None:
                missing_vals_dict[condition][tissue]['cerad']+=1
            else:
                info_dict[condition][tissue][subs]['cerad'].append(cerad)

    for rid in q_sub:
        condition,age,apoe,braak,cerad,sex,tissue = rid[1],rid[2],rid[3],rid[4],rid[5],rid[6],rid[7]
        if not condition in info_dict:
            info_dict[condition] = {}
        if not condition in missing_vals_dict:
            missing_vals_dict[condition] = {}
        if not tissue in info_dict[condition]:
            info_dict[condition][tissue] = {}
            info_dict[condition][tissue]['90+'] = {}
            info_dict[condition][tissue]['under_90'] = {}
        if not tissue in missing_vals_dict[condition]:
            missing_vals_dict[condition][tissue] = {}
            for i in sub_keys:
                info_dict[condition][tissue]['90+'][i] = []
                info_dict[condition][tissue]['under_90'][i] = []
                if not i in missing_vals_dict[condition][tissue]:
                    missing_vals_dict[condition][tissue][i] = 0
        if age == '90+':
            subs = '90+'
            if sex == None:
                missing_vals_dict[condition][tissue]['sex']+=1
            else:
                info_dict[condition][tissue][subs]['sex'].append(sex)
            info_dict[condition][tissue][subs]['age'].append(age)
            if apoe == None:
                missing_vals_dict[condition][tissue]['apoe']+=1
            else:
                info_dict[condition][tissue][subs]['apoe'].append(apoe)
            if braak == None:
                missing_vals_dict[condition][tissue]['braak']+=1
            else:
                info_dict[condition][tissue][subs]['braak'].append(braak)
            if cerad == None:
                missing_vals_dict[condition][tissue]['cerad']+=1
            else:
                info_dict[condition][tissue][subs]['cerad'].append(cerad)
        if age == None:
            missing_vals_dict[condition][tissue]['age']+=1
        if age != '90+' and age != None:
            subs = 'under_90'
            if sex == None:
                missing_vals_dict[condition][tissue]['sex']+=1
            else:
                info_dict[condition][tissue][subs]['sex'].append(sex)
            info_dict[condition][tissue][subs]['age'].append(age)
            if apoe == None:
                missing_vals_dict[condition][tissue]['apoe']+=1
            else:
                info_dict[condition][tissue][subs]['apoe'].append(apoe)
            if braak == None:
                missing_vals_dict[condition][tissue]['braak']+=1
            else:
                info_dict[condition][tissue][subs]['braak'].append(braak)
            if cerad == None:
                missing_vals_dict[condition][tissue]['cerad']+=1
            else:
                info_dict[condition][tissue][subs]['cerad'].append(cerad)

    
    for con in info_dict:
        for tissue in info_dict[con]:
            ## AGE
            over_list = info_dict[con][tissue]['90+']['age']
            under_list = info_dict[con][tissue]['under_90']['age']
            sorted_ages = sorted(under_list)
            min_age = str(min(sorted_ages))
            min_age_string = f"{min_age} - 89.9"
            max_age_string = "90+"
            
            ## SEX
            over_male,under_male,over_female,under_female = 0,0,0,0
            total_sex = len(info_dict[con][tissue]['90+']['sex']) + len(info_dict[con][tissue]['under_90']['sex'])
            for s in info_dict[con][tissue]['90+']['sex']:
                if s == 1:
                    over_male+=1
                if s == 0:
                    over_female+=1
            for se in info_dict[con][tissue]['under_90']['sex']:
                if se == 1:
                    under_male+=1
                if se == 0:
                    under_female+=1
            over_f_percent,over_m_percent,under_f_percent,under_m_percent = f"{round((over_female/total_sex) *100)}%",f"{round((over_male/total_sex) * 100)}%",f"{round((under_female/total_sex) * 100)}%",f"{round((under_male/total_sex) * 100)}%"

            ## APOE
            over_e4e4,over_e3e4,over_e3e3,over_e2e4,over_e2e3,over_e2e2 = 0,0,0,0,0,0
            under_e4e4,under_e3e4,under_e3e3,under_e2e4,under_e2e3,under_e2e2 = 0,0,0,0,0,0
            over_apoes = info_dict[con][tissue]['90+']['apoe']
            for o in over_apoes:
                if o == 22.0:
                    over_e2e2+=1
                elif o == 23.0:
                    over_e2e3+=1
                elif o == 24.0:
                    over_e2e4+=1
                elif o == 33.0:
                    over_e3e3+=1
                elif o == 34.0:
                    over_e3e4+=1
                else:
                    if o == 44.0:
                        over_e4e4+=1
            under_apoes = info_dict[con][tissue]['under_90']['apoe']
            for u in under_apoes:
                if u == 22.0:
                    under_e2e2+=1
                elif u == 23.0:
                    under_e2e3+=1
                elif u == 24.0:
                    under_e2e4+=1
                elif u == 33.0:
                    under_e3e3+=1
                elif u == 34.0:
                    under_e3e4+=1
                else:
                    if u == 44.0:
                        under_e4e4+=1
                
            ## BRAAK
            over_one,over_two,over_three = 0,0,0
            under_one,under_two,under_three = 0,0,0
            over_braaks = info_dict[con][tissue]['90+']['braak']
            for ob in over_braaks:
                if ob == 0 or ob == 1 or ob == 2:
                    over_one+=1
                elif ob == 3 or ob == 4:
                    over_two+=1
                else:
                    if ob == 5 or ob == 6:
                        over_three+=1
        
            under_braaks = info_dict[con][tissue]['under_90']['braak']
            for ub in under_braaks:
                if ub == 0 or ub == 1 or ub == 2:
                    under_one+=1
                elif ub == 3 or ub == 4:
                    under_two+=1
                else:
                    if ub == 5 or ub == 6:
                        under_three+=1

            ## CERAD
            over_c1,over_c2,under_c1,under_c2 = 0,0,0,0
            over_cerads = info_dict[con][tissue]['90+']['cerad']
            for oc in over_cerads:
                if oc == 1 or oc == 2:
                    over_c1+=1
                else:
                    if oc == 3 or oc == 4:
                        over_c2+=1
            under_cerads = info_dict[con][tissue]['under_90']['cerad']
            for uc in under_cerads:
                if uc  ==1 or uc == 2:
                    under_c1+=1
                else:
                    if uc == 3 or uc == 4:
                        under_c2+=1 

            print(f"{con}:\n {tissue}: \n {min_age_string}: \n AGE: ({len(under_list)}), \n SEX: F: {under_female} ({under_f_percent}), M {under_male} ({under_m_percent}), \n APOE: E2E2 = {under_e2e2},E2E3 = {under_e2e3}, E2E4 = {under_e2e4}, E3E3 = {under_e3e3}, E3E4 = {under_e3e4} ,E4E4 = {under_e4e4},  \n BRAAK: 0 -II = {under_one} , III-IV = {under_two}, V-VI = {under_three} \n CERAD: 1-2 = {under_c1},3-4 = {under_c2} \n *************************\n {max_age_string}: \n AGE: ({len(over_list)}), \n SEX: F {over_female} ({over_f_percent}), M {over_male} ({over_m_percent}) \n APOE: E2E2 = {over_e2e2},E2E3 = {over_e2e3}, E2E4 = {over_e2e4}, E3E3 = {over_e3e3}, E3E4 = {over_e3e4} ,E4E4 = {over_e4e4},  \n BRAAK: 0 -II = {over_one} , III-IV = {over_two}, V-VI = {over_three} \n CERAD: 1-2 = {over_c1},3-4 = {over_c2}")
            if missing_vals_dict:
                missing_ages = missing_vals_dict[con][tissue]['age']
                missing_sexes = missing_vals_dict[con][tissue]['sex']
                missing_apoes = missing_vals_dict[con][tissue]['apoe']
                missing_braaks = missing_vals_dict[con][tissue]['braak']
                missing_cerads = missing_vals_dict[con][tissue]['cerad']
                print(f"{con} in {tissue} is missing: \n Age: {missing_ages}, Sex: {missing_sexes}, Apoe: {missing_apoes}, Braak: {missing_braaks}, Cerad: {missing_cerads}")
    return

def get_boxplot_dat(db):
    outfile = open("BB_AD_CON_DAT_for_boxplot.csv",'w')
    out_write = csv.writer(outfile)
    header = ['id','type','condition','count']
    out_write.writerow(header)

    out_list = []

    tab_name = "All_Info_TopGenes_tab"
    d = sql.connect(db)
    d_con = d.cursor()

    query = f"SELECT rids,a1type,con, COUNT(*) AS total_count FROM {tab_name} GROUP BY rids,a1type,con ORDER BY rids,a1type,con"

    res = d_con.execute(query).fetchall()

    for r in res:
        out_row = [r[0],r[1],r[2],r[3]]
        out_list.append(out_row)

    for rows in out_list:
        out_write.writerow(rows)
    return


def check_tissue_matches(db,tab_3):
    ids_list,bb_list = [],[]
    first_tab,second_tab = "All_Info_TopGenes_tab","SecondSet_AllInfo_tab"
    match_dict = {}
    tissue_spec_dict = {}
    tlist = ['bulkbrain','bulkbrain','bulkbrain','blood','blood','blood']
    clist = ['Control','MCI','AD','Control','MCI','AD']
    counts = [30,27,29,89,38,34]
    for ti in range(0,len(tlist)):
        tissues = tlist[ti]
        conditions = clist[ti]
        countings = counts[ti]
        if not tissues in tissue_spec_dict.keys():
            tissue_spec_dict[tissues] = {}
        tissue_spec_dict[tissues][conditions] = countings
    tab3 = pd.read_csv(tab_3)
    ensgs = pd.unique(tab3['ensg'].tolist())
    db_conn = sql.connect(db)
    d = db_conn.cursor()
    for en in ensgs:
        match_dict[en] = {}
        match_dict[en]['gene'] = ''
        posses = tab3[tab3['ensg'] == en]['pos'].unique().tolist()   
        for p in posses:
            match_dict[en][p] = {}
            tis_set = list(set(tlist))
            cons_set = list(set(clist))
            for tis in tis_set:
                match_dict[en][p][tis] = {}
                for cons in cons_set:
                    match_dict[en][p][tis][cons] = ''
                    if tis == 'bulkbrain':
                        if cons == 'AD' or cons == 'Control':
                            bb_first = f"""SELECT gene_name,con,COUNT(DISTINCT rids) FROM {first_tab} WHERE ref_base = ? AND poss = ? GROUP BY con"""
                            bb = d.execute(bb_first,(en,p)).fetchall()
                            gene_name = bb[0]
                            if not gene_name == match_dict[en]['gene']:
                                match_dict[en]['gene'] = gene_name
                            bb_list.append(bb)
                        if cons == 'MCI':
                            bb_second = f"""SELECT ref_id,condition,COUNT(DISTINCT rid) FROM {second_tab} WHERE ref_id = ? AND pos = ? AND condition = ? AND tissue = ? GROUP BY condition"""
                            bb2 = d.execute(bb_second,(en,p,cons,tis)).fetchall()
                            bb_list.append(bb2)
    print(bb_list)

    return


def build_profiles(db,in_file,samp_json,profiles_json,profile_key_json):
    first_tab,second_tab = "All_Info_TopGenes_tab","SecondSet_AllInfo_tab"
    
    db_conn = sql.connect(db)
    d = db_conn.cursor()
    lists  = pd.read_csv(in_file)
    genes = lists[['gene','type','position','ensg']].drop_duplicates()
    sets = [tuple(row)for row in genes.to_numpy()]
    ct = 0
    
    make_profiles_dict,profile_key_dict = {},{}
    samp_dict = {}
    for i in ['Control','AD','MCI']:
        samp_dict[i] = {}
        for j in ['bulkbrain','monocyte']:
            samp_dict[i][j] = {}
    for se in sets:
        ct+=1
        c = str(ct)
        gene,types,poss,ensg = se[0],se[1],se[2],se[3]
        mod = f"mod_{c}"
        profile_key_dict[mod] = se
        rids1 = f"SELECT rids,con,age_death,msex,apoe_genotype,braaksc,ceradsc,ref,alt1 FROM {first_tab} WHERE gene_name = ? AND poss = ?"
        rids2 = f"SELECT rid,condition,age,sex,apoe,braak,cerad,tissue,ref,alt1 FROM {second_tab} WHERE ref_id = ? AND pos = ?"
        r1 = d.execute(rids1,(gene,poss)).fetchall()
        r2 = d.execute(rids2,(ensg,poss)).fetchall()

        if r1:
            for ones in r1:
                samp,condition,age,sex,apoe,braak,cerad,tissue,ref,alt1 = ones[0],ones[1],ones[2],ones[3],ones[4],ones[5],ones[6],'bulkbrain',ones[7],ones[8]
                if not samp in samp_dict[condition][tissue]:
                    samp_dict[condition][tissue][samp] = {}
                    samp_dict[condition][tissue][samp]['info'] = (age,sex,apoe,braak,cerad,tissue)
                    samp_dict[condition][tissue][samp]['mods'] = []
                samp_dict[condition][tissue][samp]['mods'].append(mod)
        if r2:
            for twos in r2:
                samp,condition,age,sex,apoe,braak,cerad,tissue,ref,alt1 = twos[0],twos[1],twos[2],twos[3],twos[4],twos[5],twos[6],twos[7],twos[8],twos[9]
                if not samp in samp_dict[condition][tissue]:
                    samp_dict[condition][tissue][samp] = {}
                    samp_dict[condition][tissue][samp]['info'] = (age,sex,apoe,braak,cerad,tissue)
                    samp_dict[condition][tissue][samp]['mods'] = []
                samp_dict[condition][tissue][samp]['mods'].append(mod)
    
    with open(profile_key_json,'w') as key_json:
        json.dump(profile_key_dict,key_json)

    for con in samp_dict.keys():
        if not con in make_profiles_dict.keys():
            make_profiles_dict[con] = {}
        for tiss in samp_dict[con].keys():
            if not tiss in make_profiles_dict[con].keys():
                make_profiles_dict[con][tiss] = {}
            for samples in samp_dict[con][tiss].keys():
                mods = tuple(samp_dict[con][tiss][samples]['mods'])
                infos = samp_dict[con][tiss][samples]['info']
                if not str(mods) in make_profiles_dict[con][tiss].keys():
                    make_profiles_dict[con][tiss][str(mods)] = {}
                    infos_list = ['age','sex','apoe','braak','cerad']
                    for s in infos_list:
                        if not s in  make_profiles_dict[con][tiss][str(mods)]:
                            make_profiles_dict[con][tiss][str(mods)][s] = {}
                ag = infos[0]
                if type(ag) == int or type(ag) == float:
                    if ag <= 89.9:
                        ag = '<= 89.9'
                    else: 
                        ag = str(infos[0]) 
                subkeys = [ag,str(infos[1]),str(infos[2]),str(infos[3]),str(infos[4])]
                for su in range(0,len(subkeys)):
                    supkey = infos_list[su]
                    if not subkeys[su] in make_profiles_dict[con][tiss][str(mods)].keys():
                        make_profiles_dict[con][tiss][str(mods)][supkey] = {}
                        make_profiles_dict[con][tiss][str(mods)][supkey][subkeys[su]] = 0
                    make_profiles_dict[con][tiss][str(mods)][supkey][subkeys[su]]+=1

    with open(profiles_json,'w') as profs_json:
        json.dump(make_profiles_dict,profs_json)

                
    for key in make_profiles_dict:
        for val in make_profiles_dict[key]:
            for modis in make_profiles_dict[key][val]:
                print(f"Profiles build for: {key}: {val}: {modis}          {make_profiles_dict[key][val][modis]}")

    return 

def top_vars_in_monocytes(db,varsin):
    second_tab = "SecondSet_AllInfo_tab"
    seconds = {}
    final_results = []
    lists  = pd.read_csv(varsin)
    genes = lists[['ensg','pos']].drop_duplicates()
    sets = [tuple(row)for row in genes.to_numpy()]
    with futures.ProcessPoolExecutor(max_workers=32) as mst:
        print(f'************Getting counts from {second_tab}',datetime.now())
        wait_for = [mst.submit(GetCounts.get_counts,db,sets[se],se,second_tab) for se in range(0,len(sets)) ]
        for fu in futures.as_completed(wait_for):
            current = fu.result()
            final_results.append(current)
    outs = []
    for fi in final_results:
        if fi:
            mods,tiss,con,ensg,pos,ct = '','','','','',''
        
            for i in fi:
                if type(i) == str:
                    mods = i
                if type(i) == tuple:
                    tiss,con,ensg,pos,ct = i[0],i[1],i[2],i[3],i[4]
                
                infos_list = [mods,tiss,con,ensg,pos,ct]
                outs.append(infos_list)
    outfile = open('second_tab_infos.csv', 'w')
    out_write = csv.writer(outfile)
    for infos in outs:
        out_write.writerow(infos)

    print("WE ARE DONE WITH THESE ANALYSES")
        
    return

if __name__=='__main__':
   dbp = sys.argv[1]
   mono_vcf = sys.argv[2]
   bb_vcf =  sys.argv[3]
   clin =  sys.argv[4]
   out = sys.argv[5]
   bb_firstset = sys.argv[6]
   mod_dump = sys.argv[7]
   inpath_list = [mono_vcf,bb_vcf]
   inpath_tissues = ['monocyte','bulkbrain']
   in_mods = sys.argv[8]
   out_prof = sys.argv[9]
   j1 = sys.argv[10]
   j2 = sys.argv[11]
   j3 = sys.argv[12]
   topvars = sys.argv[13]
   #mod_dict = make_modid_dict(bb_firstset,mod_dump)
   #secondset_tab = make_samp_tabs(dbp,inpath_list,inpath_tissues,clin,out,mod_dict)
   #results = all_tissues_conditions(dbp,out,mod_dict)
   #sizes = get_sample_sizes(dbp)
   #volcano_dat = run_analysis_for_volcano(out)
   #samples_info = get_samp_info(dbp,clin)
   #boxy = get_boxplot_dat(dbp)
   tissue_check = check_tissue_matches(dbp,bb_firstset)
   #profiles = build_profiles(dbp,in_mods,j2,j1,j3)
   #monocytes = top_vars_in_monocytes(dbp,topvars)
