import csv
import os
import sys
import sqlite3 as sql
import pandas as pd
from scipy.stats import mannwhitneyu
from collections import Counter
from scipy.stats import ttest_ind
from position_groups import PosGroups

def get_sample_info(db):
    keys_list,vals_list = ['braak','cerad','apoe','sex'],[['V-VI','III-IV','I-II','None'],['1-2','3-4'],['E4E4','E3E4','E3E3','E2E4','E2E3','E2E2','None'],['m','f']]
    
    val2_list = [['6.0','5.0','4.0','3.0','2.0','1.0','0.0'],['1.0','2.0','3.0','4.0'],['44.0','34.0','33.0','24.0','23.0','22.0','None'],['1.0','0.0']]
    
    braak_idx,cerad_idx,apoe_idx,sex_idx = keys_list.index('braak'),keys_list.index('cerad'),keys_list.index('apoe'),keys_list.index('sex')

    samples = {}
    samples['AD'] = {}
    samples['Control'] = {}
    for key in samples.keys():
        samples[key]['90+']={}
        samples[key]['<90']={}
        for k1 in samples[key].keys():
            samples[key][k1] = {}
            for l in range(0,len(keys_list)):
                k2 = keys_list[l]
                v_list = vals_list[l]
                samples[key][k1][k2] = {}
                for va in v_list:
                    samples[key][k1][k2][va]=0
    conn = sql.connect(db)
    c = conn.cursor()
    table = 'All_Info_TopGenes_tab'
    
    query = c.execute(f"SELECT DISTINCT rids, apoe_genotype, braaksc,ceradsc, msex, age_death, con FROM {table}").fetchall()
    
    for qu in query:
        apoe,braak,cerad,sex,age,con = str(qu[1]),str(qu[2]),str(qu[3]),str(qu[4]),str(qu[5]),str(qu[6])
        print(con,age)
        braaks,cerads,apoes,sexes = vals_list[braak_idx],vals_list[cerad_idx],vals_list[apoe_idx],vals_list[sex_idx]
        braak_key,cerad_key,apoe_key,sex_key = '','','','',
        if braak in val2_list[braak_idx][0:2]:
            braak_key = braaks[0]
        if braak in val2_list[braak_idx][2:4]:
            braak_key = braaks[1]
        if braak in val2_list[braak_idx][4:6]:
            braak_key = braaks[2]
        if braak in val2_list[braak_idx][-1]:
            braak_key = braaks[3]
        
        if cerad in val2_list[cerad_idx][0:2]:
            cerad_key = cerads[0]
        if cerad in val2_list[cerad_idx][2:4]:
            cerad_key = cerads[1]

        if not apoe == 'None':
            apo = val2_list[apoe_idx].index(apoe)
            apoe_key = apoes[apo]
        else:
            apoe_key = 'None'

        se = val2_list[sex_idx].index(sex)
        sex_key = sexes[se]
        
         
        subs_list =  [braak_key,cerad_key,apoe_key,sex_key]
        if age == '90+':
            for r in range(0,len(keys_list)):
                k = keys_list[r]
                rv = subs_list[r]
                samples[con]['90+'][k][rv] +=1
        if age != '90+':
            for r in range(0,len(keys_list)):
                k = keys_list[r]
                rv = subs_list[r]
                samples[con]['<90'][k][rv] +=1


    for conditions in samples:
        print(f"Condition: {conditions}")
        for ages in samples[conditions].keys():
            print(f"\tAge: {ages}")
            for var in samples[conditions][ages].keys():
                print(f"\t\t{var}")
                for subvar in samples[conditions][ages][var].keys():
                    ct = samples[conditions][ages][var][subvar]
                    print(f"\t\t\t{subvar} = {ct}\n")
    return



def check_distinct_by_edit_type(db):
    con = sql.connect(db)
    c = con.cursor()
    table = 'All_Info_TopGenes_tab'
    query = c.execute(f"SELECT a1type, COUNT(DISTINCT poss) AS distinct_edit_type_count FROM {table} GROUP BY a1type").fetchall()
    for q in query:
        print(q)
    return

def table2_info(db):
    con = sql.connect(db)
    c = con.cursor()
    table = 'All_Info_TopGenes_tab'
    unique_check = {}
    uniques = {}
    commons = {}

    ##For unique vs common edits
    query = c.execute(f"SELECT DISTINCT refid,poss,COUNT(rids) AS cpx FROM {table} GROUP BY refid,poss").fetchall()
    for q in query:
        refid,poss,ct = q[0],str(q[1]),q[2]
        if not refid in unique_check.keys():
            unique_check[refid] = {}
        if not poss in unique_check[refid].keys():
            unique_check[refid][poss]=0
        unique_check[refid][poss]+=ct
    for refs in unique_check.keys():
        for posi in unique_check[refs].keys():
            cts = unique_check[refs][posi]
            if cts == 1:
                if not refs in uniques.keys():
                    uniques[refs] = []
                uniques[refs].append(posi)
            else:
                if not refs in commons.keys():
                    commons[refs] = []
                commons[refs].append(posi)
    ##Get deets for uniques
    #uniques_keys,uniques_values = [],[]
    #for ke in uniques.keys():
    #    uniques_keys.append(ke)
    #for li in uniques.values():
    #    for i in li:
    #        uniques_values.append(i)
    #keys_ph = ','.join('?' * len(uniques_keys))
    #val_ph = ','.join('?' * len(uniques_values))
    #query2 = f"""SELECT DISTINCT edit_dist,a1type,con, COUNT(poss) FROM {table} WHERE refid IN ({keys_ph}) and poss IN ({val_ph}) GROUP BY edit_dist,a1type,con"""
    #params = tuple(uniques_keys) + tuple(uniques_values)
    #quer2 = c.execute(query2,params).fetchall()
    #uniques_outs = {}
    #for q2 in quer2:
    #    dist,types,cons,cts = q2[0],q2[1],q2[2],q2[3]
    #    if "AD" in dist:
    #        dist = "AD_enriched"
    #    if "Control" in dist:
    #        dist = "Control_enriched"
    #    if not dist in uniques_outs.keys():
    #        uniques_outs[dist] = {}
    #    if not types in uniques_outs[dist].keys():
    #        uniques_outs[dist][types] = {}
    #    if not cons in  uniques_outs[dist][types].keys():
    #        uniques_outs[dist][types][cons]=0
    #    uniques_outs[dist][types][cons]+=cts
    #for d in uniques_outs.keys():
    #    for t in uniques_outs[d].keys():
    #        for c in uniques_outs[d][t].keys():
    #            un_cts = uniques_outs[d][t][c]
    #            print(f"UNIQUE: {d}\t{t}\t{c}\t{str(un_cts)}")
    ##Get deets for commons
    commons_keys,commons_values = [],[]
    for ke in commons.keys():
        commons_keys.append(ke)
    for li in commons.values():
        for i in li:
            commons_values.append(i)
    keys_ph = ','.join('?' * len(commons_keys))
    val_ph = ','.join('?' * len(commons_values))
    query3 = f"""SELECT DISTINCT edit_dist,a1type,con, COUNT(DISTINCT poss) FROM {table} WHERE refid IN ({keys_ph}) and poss IN ({val_ph}) GROUP BY edit_dist,a1type,con"""
    params = tuple(commons_keys) + tuple(commons_values)
    quer3 = c.execute(query3,params).fetchall()
    commons_outs = {}
    for q3 in quer3:
        dist,types,cons,cts = q3[0],q3[1],q3[2],q3[3]
        if "AD" in dist:
            dist = "AD_enriched"
        if "Control" in dist:
            dist = "Control_enriched"
        if not dist in commons_outs.keys():
            commons_outs[dist] = {}
        if not types in commons_outs[dist].keys():
            commons_outs[dist][types] = {}
        if not cons in  commons_outs[dist][types].keys():
            commons_outs[dist][types][cons]=0
        commons_outs[dist][types][cons]+=cts
    for d in commons_outs.keys():
        for t in commons_outs[d].keys():
            for c in commons_outs[d][t].keys():
                com_cts = commons_outs[d][t][c]
                print(f"COMMON: {d}\t{t}\t{c}\t{str(com_cts)}")
    return

def global_local_edit_types(db):
    con = sql.connect(db)
    c = con.cursor()
    table = 'All_Info_TopGenes_tab'
    global_query = c.execute(f"SELECT DISTINCT a1type,con,rids,SUM(a1prop) FROM {table} GROUP BY a1type,con,rids").fetchall()
    global_dict = {}
    for gq in global_query:
        etype,cond,sums = gq[0],gq[1],gq[3]
        if not etype in global_dict.keys():
            global_dict[etype] = {}
        if not cond in global_dict[etype].keys():
            global_dict[etype][cond] = []
        global_dict[etype][cond].append(sums)
    for edit_type in global_dict.keys():
        ads,cons = [],[]
        ad,co = 29,30
        for cono in global_dict[edit_type].keys():
            pros = global_dict[edit_type][cono]
            le = len(pros)
            if cono == "AD":
                for pro in pros:
                    ads.append(pro)
                diff = ad-le
                ads.extend([0]*diff)
            if cono == "Control":
                for pro in pros:
                    cons.append(pro)
                diff = co-le
                cons.extend([0]*diff)
        sig = ''
        stat, p_value = mannwhitneyu(ads, cons, alternative='two-sided')
        if float(p_value) <= float(0.05):
            sig = "is"
        else:
            sig = "is not"

        #print(f"{edit_type}: global diff p-val = {float(p_value)} {sig} significant at threshold 0.05")    
    
    local_dict = {}
    g_list = []
    cols = ['refid','type','sig_lvl','num_types','ad_ct','con_ct']
    local_df = pd.DataFrame(columns=cols)
    #local_file = open("local_level_highly_sig_type_diffs_per_gene.csv",'w')
    #local_write = csv.writer(local_file)
    local_query = c.execute(f"SELECT DISTINCT refid,a1type,p_value,con, COUNT(poss) FROM {table} GROUP BY refid,a1type,con").fetchall()
    for loc in local_query:
        g,t,p,c,ct = loc[0],loc[1],loc[2],loc[3],loc[4]
        g = '_'.join(g.split('_')[:-1])
        
        if p:
            p = float(p)
            sp = ''
            if p < float(0.001):
                sp = '***'
            elif p < float(0.01):
                sp = '**'
            else:
                if p < float(0.05):
                    sp = '*'
        
            out_tupe = (g,sp,c,ct)
            
            g_list.append(g)
            if not t in local_dict.keys():
                local_dict[t] = []
            if not out_tupe in local_dict[t]:
                local_dict[t].append(out_tupe)

    tots = Counter(g_list)
    to = dict(tots)
    outs = []
    ##ADDS NUMBER (how many other edit-types signif for that gene?) TO END OF GENE ID
    for types in local_dict.keys():
        new_dict = {}
        for tupes in local_dict[types]:
            genes,sigs,cons,cts = tupes[0],tupes[1],tupes[2],tupes[3]
            num = 0
            if genes in to.keys():
                num = to[genes]
            if not genes in new_dict.keys():
                new_dict[genes] = {}
            if not sigs in new_dict[genes].keys():
                new_dict[genes][sigs] = {}
                new_dict[genes][sigs]['Control'] = None
                new_dict[genes][sigs]['AD'] = None
                new_dict[genes][sigs]['num'] = num
            new_dict[genes][sigs][cons] = cts
            
           
        for ge in new_dict.keys():
            for si in new_dict[ge].keys():
                ad_ct = new_dict[ge][si]['AD']
                con_ct = new_dict[ge][si]['Control']
                nu  = new_dict[ge][si]['num']
                outrow = [ge,types,si,nu,ad_ct,con_ct]
                outs.append(outrow)
        
    for ro in outs:
        local_df.loc[len(local_df)] = ro

    filtered = local_df[(local_df['sig_lvl'] == "***") & ((local_df['ad_ct'] >= 15) | (local_df['con_ct'] >= 15))]
    filtered.to_csv("local_level_highly_sig_type_diffs_per_gene.csv",index = False)
    
    
    local_df.to_csv("local_level_all_sig_type_diffs_per_gene.csv",index=False)
    
    
    unique_local_genes = filtered['refid'].unique().tolist()
    print(unique_local_genes,len(unique_local_genes))
    
    return

def local_set_deets(db,flt):
    out_df = pd.DataFrame(columns=['refid','edit_type','pos','condition','samp','a1ad','dp','apoe_genotype','braak','cerad','hgvs','ref','alt1'])
    i = 0
    filt = pd.read_csv(flt)
    filt['n_ad_samps'] = None
    filt['n_con_samps'] = None
    id_type_list = []
    con = sql.connect(db)
    c = con.cursor()
    table = 'All_Info_TopGenes_tab'
    unique_local_genes = filt['refid'].unique().tolist()
    dist_samps = c.execute(f"SELECT DISTINCT rids,con FROM {table}").fetchall()
    #columns = [info[1] for info in c.fetchall()]
    all_samps = []
    for dist in dist_samps:
        s,condi = dist[0],dist[1]
        out_s = f"{s}_{condi}"
        all_samps.append(out_s)
    for re in unique_local_genes:
        temp_df = pd.DataFrame(columns=['refid','edit_type','pos','condition','samp','a1ad','dp','apoe_genotype','braak','cerad','edit_dist','ref','alt1'])
        df_list = []
        types = filt[filt['refid'] == re][['type']]
        for ty in types['type']:
            q_state = f"SELECT poss,con,rids,a1ad,dp,apoe_genotype,braaksc,ceradsc,edit_dist,ref,alt1 FROM {table} WHERE ref_base = ? AND a1type = ?"
            qs = c.execute(q_state,(re,ty)).fetchall()
            for q in qs:
                poss,con,rids,a1ad,dp,apoe_genotype,braaksc,ceradsc,edist,ref,alt1 = q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10]
                out_row = [re,ty,poss,con,rids,a1ad,dp,apoe_genotype,braaksc,ceradsc,edist,ref,alt1]
                df_list.append(out_row)
        for lists in df_list:
            temp_df.loc[len(temp_df)] = lists
        
        if i == 0:
            out_df = pd.concat([out_df, temp_df], ignore_index=True)
        else:
            out_df = pd.concat([out_df, temp_df], ignore_index=True)
        i+=1

        ## Count unique samples per edit type per gene
       
        tys = temp_df['edit_type'].unique().tolist()
        for typ in tys:

            refs = temp_df.loc[temp_df['edit_type'] == typ,'ref'].tolist()
            alt1s = temp_df.loc[temp_df['edit_type'] == typ, 'alt1'].tolist()
            print(typ,refs,alt1s)    
            #valu_list = []
            #pos_cols = ['sample']
            #out_pos = f'{re}_{typ}_positions_mat.csv'
            #out_dpth = f'{re}_{typ}_depth_mat.csv'
            #ads = temp_df[(temp_df['condition'] == 'AD') & (temp_df['edit_type'] == typ)]
            #controls = temp_df[(temp_df['condition'] == 'Control') & (temp_df['edit_type'] == typ)]
            #n_ads,n_controls = 0,0
            #n_ads_for_gene,n_controls_for_gene = [],[]
            
            #if not ads.empty:
            #    n_ads_for_gene = ads['samp'].unique().tolist()
            #    n_ads = len(n_ads_for_gene)
            #if not controls.empty:
            #    n_controls_for_gene = controls['samp'].unique().tolist()
            #    n_controls = len(n_controls_for_gene)
            
            #mask = (filt['refid'] == re) & (filt['type'] == typ)
            #filt.loc[mask, 'n_con_samps'] = n_controls
            #filt.loc[mask, 'n_ad_samps'] = n_ads
             
            #positions = temp_df['pos'].unique().tolist()
            #positions_list = []
            #for po in positions:
            #    po = f"{po}_{typ}"
            #    pos_cols.append(po)
            #pos_df = pd.DataFrame(columns = pos_cols)
            #pos_dict_list = []
            #valu_list = []
            #depth_df = pd.DataFrame(columns = ['refid','edit_type','samp','con','pos','a1ad','dp','edit_dist','apoe_genotype','braak','cerad','sex'])
            #pos_df_list = []
            #for sam in all_samps:
            #    posit_list = []
            #    pos_dict = {}
            #    if not sam in pos_dict.keys():
            #        pos_dict[sam] = {}
            #    sa = sam.split('_')[0]
                
            #    for posi in positions:
            #        posit = f"{posi}_{typ}"
            #        pos_cols.append(posit)
            #        posit_list.append(posit)
            #        if not posit in pos_dict[sam].keys():
            #            pos_dict[sam][posit] = ''
            #        sm_pos = temp_df[(temp_df['samp'] == sa) & (temp_df['pos'].astype(int) == int(posi)) & (temp_df['edit_type'] == typ) & (temp_df['refid']== re)]
        
            #        if not sm_pos.empty:
            #            pos_dict[sam][posit] = 'yes'
            #            valu_list.append(sm_pos[['refid','edit_type','samp','condition','pos','a1ad','dp','edit_dist']].values.tolist())
            #        else:
            #            pos_dict[sam][posit] = 'no'
            #    for sss in pos_dict.keys():
            #        s_list = [sss]
            #        for poi in pos_dict[sss].keys():
            #            decision = pos_dict[sss][poi]
            #            s_list.append(decision)
            #        pos_df.loc[len(pos_df)] = s_list
            
            
            #pos_df.to_csv(out_pos, index=False)
            #for vals in valu_list:
            #    row = vals[0]
            #    the_samp = row[2]
            #    get_extra = f"SELECT DISTINCT apoe_genotype,braaksc,ceradsc,msex FROM {table} where rids = ?"
            #    gets = c.execute(get_extra,(the_samp,)).fetchall()
            #    for get in gets:
            #        apoe_geno,bra,cer,ms = get[0],get[1],get[2],get[3]
            #        for inf in [apoe_geno,bra,cer,ms]:
            #            row.append(inf)
            #    depth_df.loc[len(depth_df)] = row

            #ad_depths,ad_perc = [],[]
            #control_depths,control_perc = [],[]

            #check_ad,check_con = {},{}

            #add = depth_df[depth_df['con'] == "AD"][['a1ad','dp','samp','pos']].values.tolist()
           # for al in add:
           #     ade,dps,ats,ap = al[0],al[1],al[2],al[3]
                
           #     if not ats in check_ad.keys():
           #         check_ad[ats] = {}
           #     if not ap in check_ad[ats].keys():
           #         check_ad[ats][ap] = None
           #     check_ad[ats][ap] = (ade,dps)
           # for all_ad in check_ad.keys():
           #     all_ads,all_dep = 0,0
           #     for posing in check_ad[all_ad].keys():
           #         out_tupe = check_ad[all_ad][posing]
           #         alla,alld = out_tupe[0],out_tupe[1]
           #         all_dep+=alld
           #         all_ads+=alla
           #     ad_depths.append(all_dep)
           #     a_p = int((all_ads/all_dep)*100)
           #     ad_perc.append(a_p)
            

           # cond = depth_df[depth_df['con'] == "Control"][['a1ad','dp','samp','pos']].values.tolist()
           # for cl in cond:
           #     cde,dps,cts,cp = cl[0],cl[1],cl[2],cl[3]

           #     if not cts in check_con.keys():
           #         check_con[cts] = {}
           #     if not cp in check_con[cts].keys():
           #         check_con[cts][cp] = None
           #     check_con[cts][cp] = (cde,dps)
           # for cll_control in check_con.keys():
           #     cll_ads,cll_dep = 0,0
           #     for posingc in check_con[cll_control].keys():
           #         out_tupe = check_con[cll_control][posingc]
           #         clla,clld = out_tupe[0],out_tupe[1]
           #         cll_dep+=clld
           #         cll_ads+=clla
           #     control_depths.append(cll_dep)
           #     c_p = int((cll_ads/cll_dep)*100)
           #     control_perc.append(c_p)
           # if len(ad_depths) > 1 and len(control_depths) > 1:
           #     t_stat_depths, p_value_depths = ttest_ind(ad_depths, control_depths, equal_var=False)
           #     t_stat_percs, p_value_percs = ttest_ind(ad_perc, control_perc, equal_var=False)
           # print(f"Difference in read depths between AD and Control for edit-type {typ} in geneid {re} returned a p-val of {p_value_depths}")
           # print(f"Difference in percent gene edited between AD and Control for edit-type {typ} in geneid {re} returned a p-val of {p_value_percs}")
           # depth_df.to_csv(out_dpth,index=False)
    
            
    return 

def per_gene_investigate(db,flt):
    con = sql.connect(db)
    co = con.cursor()
    table = 'All_Info_TopGenes_tab'
    f = pd.read_csv(flt)
    refid_list = f['refid'].unique().tolist()
    ad,control = 29,30
    final_outs = "Table3_results.csv"
    final_open = open(final_outs,'w')
    final_write = csv.writer(final_open)
    header_final = ['ensg','transcriptome_pos','ad_ct','con_ct','edit_type','global_dist','local_dist','hgvs','prop_diff_pval','higher_avg_prop']
    final_write.writerow(header_final)
    #for i in ["specific","double","differential"]
    for refs in refid_list:
        dict1 = {}
        hgvs_list = []
        st = f"{refs}_hgvs.csv"
        hg_out = open(st, 'w')
        headers = ["pos","hgvs"]
        hg_out_write = csv.writer(hg_out)
        hg_out_write.writerow(headers)
        types = f[f['refid'] == refs]['type'].unique()
        typing = types.tolist()
        hgvs_dict = {}
        diff_dict = {}
        for ty in typing:
            dict1[ty] = {}
            for i in ["specific","double","differential","low-specific"]:
                dict1[ty][i] = {}
            query1 = f"SELECT con, COUNT(DISTINCT rids) FROM {table} WHERE ref_base = ? AND a1type = ? GROUP BY con"
            que = co.execute(query1, (refs,ty)).fetchall()
            ##Cndition specific
            if len(que) == 1:
                con,ct = que[0][0],que[0][1]
                if con == "AD"and ct == ad:
                    query2 = f"SELECT DISTINCT con, rids, hgvs, poss FROM {table} WHERE ref_base = ? and a1type = ?"
                    q2 = co.execute(query2,(refs,ty)).fetchall()
                    for qs in q2:
                        condition,rids,hgvs,pos = qs[0],qs[1],qs[2],qs[3]
                        hgvs_tupe = (pos,hgvs)
                        if not pos in hgvs_dict.keys():
                            hgvs_dict[str(pos)] = hgvs
                        if not hgvs_tupe in hgvs_list:
                            hgvs_list.append(hgvs_tupe)
                        if not condition in dict1[ty]["specific"].keys():
                            dict1[ty]["specific"][condition] = {}
                        if not rids in dict1[ty]["specific"][condition].keys():
                            dict1[ty]["specific"][condition][rids] = []
                        dict1[ty]["specific"][condition][rids].append(pos)
                    
                elif con == "Control" and ct == control:
                    query2 = f"SELECT DISTINCT con, rids, hgvs, poss FROM {table} WHERE ref_base = ? and a1type = ?"
                    q2 = co.execute(query2,(refs,ty)).fetchall()
                    for qs in q2:
                        condition,rids,hgvs,pos = qs[0],qs[1],qs[2],qs[3]
                        hgvs_tupe = (pos,hgvs)
                        if not pos in hgvs_dict.keys():
                            hgvs_dict[str(pos)] = hgvs
                        if not hgvs_tupe in hgvs_list:
                            hgvs_list.append(hgvs_tupe)
                        if not condition in dict1[ty]["specific"].keys():
                            dict1[ty]["specific"][condition] = {}
                        if not rids in dict1[ty]["specific"][condition].keys():
                            dict1[ty]["specific"][condition][rids] = []
                        dict1[ty]["specific"][condition][rids].append(pos)
                
                else:
                    query2 = f"SELECT DISTINCT con, rids, hgvs, poss FROM {table} WHERE ref_base = ? and a1type = ?"
                    q2 = co.execute(query2,(refs,ty)).fetchall()
                    for qs in q2:
                        condition,rids,hgvs,pos = qs[0],qs[1],qs[2],qs[3]
                        hgvs_tupe = (pos,hgvs)
                        if not pos in hgvs_dict.keys():
                            hgvs_dict[str(pos)] = hgvs
                        if not hgvs_tupe in hgvs_list:
                            hgvs_list.append(hgvs_tupe)
                        if not condition in dict1[ty]["low-specific"].keys():
                            dict1[ty]["low-specific"][condition] = {}
                        if not rids in dict1[ty]["low-specific"][condition].keys():
                            dict1[ty]["low-specific"][condition][rids] = []
                        dict1[ty]["low-specific"][condition][rids].append(pos)
            else:
                con1,con2 = que[0],que[1]
                condi1,condi2 = con1[0],con2[0]
                ct1,ct2 = con1[1],con2[1]
                if ct1 > ct2:
                    if ct1 >= 1.5 * ct2:
                        query2 = f"SELECT DISTINCT con, rids, hgvs, poss FROM {table} WHERE ref_base = ? and a1type = ?"
                        q2 = co.execute(query2,(refs,ty)).fetchall()
                        for qs in q2:
                            condition,rids,hgvs,pos = qs[0],qs[1],qs[2],qs[3]
                            hgvs_tupe = (pos,hgvs)
                            if not pos in hgvs_dict.keys():
                                hgvs_dict[str(pos)] = hgvs
                            if not hgvs_tupe in hgvs_list:
                                hgvs_list.append(hgvs_tupe)
                            if not condition in dict1[ty]["double"].keys():
                                dict1[ty]["double"][condition] = {}
                            if not rids in dict1[ty]["double"][condition].keys():
                                dict1[ty]["double"][condition][rids] = []
                            dict1[ty]["double"][condition][rids].append(pos)
                    else:

                        
                        query2 = f"SELECT DISTINCT con, rids, hgvs, poss,a1ad,dp FROM {table} WHERE ref_base = ? and a1type = ?"
                        q2 = co.execute(query2,(refs,ty)).fetchall()
                        for qs in q2:
                            con_list,ad_list = [],[]
                            condition,rids,hgvs,pos,adepth,dp = qs[0],qs[1],qs[2],qs[3],qs[4],qs[5]
                            prop = adepth/dp
                            if condition == "AD":
                                ad_list.append(prop)
                            if condition == "Control":
                                con_list.append(prop)
                            hgvs_tupe = (pos,hgvs)
                            if not pos in hgvs_dict.keys():
                                hgvs_dict[str(pos)] = hgvs
                            if not hgvs_tupe in hgvs_list:
                                hgvs_list.append(hgvs_tupe)
                            if not condition in dict1[ty]["differential"].keys():
                                dict1[ty]["differential"][condition] = {}
                            if not rids in dict1[ty]["differential"][condition].keys():
                                dict1[ty]["differential"][condition][rids] = []
                            dict1[ty]["differential"][condition][rids].append(pos)
                            
                if ct2 > ct1:
                    if ct2 >= 1.5 * ct1:
                        query2 = f"SELECT DISTINCT con, rids, hgvs, poss FROM {table} WHERE ref_base = ? and a1type = ?"
                        q2 = co.execute(query2,(refs,ty)).fetchall()
                        for qs in q2:
                            condition,rids,hgvs,pos = qs[0],qs[1],qs[2],qs[3]
                            hgvs_tupe = (pos,hgvs)
                            if not pos in hgvs_dict.keys():
                                hgvs_dict[str(pos)] = hgvs
                            if not hgvs_tupe in hgvs_list:
                                hgvs_list.append(hgvs_tupe)
                            if not condition in dict1[ty]["double"].keys():
                                dict1[ty]["double"][condition] = {}
                            if not rids in dict1[ty]["double"][condition].keys():
                                dict1[ty]["double"][condition][rids] = []
                            dict1[ty]["double"][condition][rids].append(pos)
                    else:
                        
                        query2 = f"SELECT DISTINCT con, rids, hgvs, poss,a1ad,dp FROM {table} WHERE ref_base = ? and a1type = ?"
                        q2 = co.execute(query2,(refs,ty)).fetchall()
                        ad_list,con_list = [],[]
                        
                        for qs in q2:
                            condition,rids,hgvs,pos,adepth,dp = qs[0],qs[1],qs[2],qs[3],qs[4],qs[5]
                            prop = adepth/dp
                            if condition == "AD":
                                ad_list.append(prop)
                            if condition == "Control":
                                con_list.append(prop)
                            hgvs_tupe = (pos,hgvs)
                            if not pos in hgvs_dict.keys():
                                hgvs_dict[str(pos)] = hgvs
                            if not hgvs_tupe in hgvs_list:
                                hgvs_list.append(hgvs_tupe)
                            if not condition in dict1[ty]["differential"].keys():
                                dict1[ty]["differential"][condition] = {}
                            if not rids in dict1[ty]["differential"][condition].keys():
                                dict1[ty]["differential"][condition][rids] = []
                            dict1[ty]["differential"][condition][rids].append(pos)
            
        
        
        ## edit-type
        for di in dict1.keys():
            ## local-dist
            for dik in dict1[di].keys():
                if dict1[di].get(dik):
                    full_dict = {}
                    full_dict[refs] = {}
                    full_dict[refs][di] = {}
                    full_dict[refs][di][dik] = {}
                    full_dict[refs][di][dik]['pos'] = {}
                    pos_list,ad_list,con_list = [],[],[]
                    for condition in dict1[di][dik].keys():
                        for ids in dict1[di][dik][condition].keys():
                            posis = dict1[di][dik][condition][ids]
                            for pos in posis:
                                pi  = str(pos)
                                if not pi in full_dict[refs][di][dik]['pos'].keys():
                                    full_dict[refs][di][dik]['pos'][pi] = {}
                                    full_dict[refs][di][dik]['pos'][pi]['ad_ct'] = 'na'
                                    full_dict[refs][di][dik]['pos'][pi]['con_ct'] = 'na'
                                    
                                    hg = hgvs_dict[pi]
                                    full_dict[refs][di][dik]['pos'][pi]['hgvs'] = hg
                                    

                                pos_tupe = (pi,condition)
                                if not pos_tupe in pos_list:
                                    
                                    condict = {}
                                    pos_list.append(pos_tupe)
                                    query1 = f"SELECT con, COUNT(DISTINCT rids) FROM {table} WHERE ref_base = ? AND a1type = ? and poss = ? GROUP BY con"
                                    que1 = co.execute(query1, (refs,di,pos)).fetchall()
                                    for q1 in que1:
                                        condit,count = q1[0],q1[1]
                                        condict[condit] = str(count) 
                                    if full_dict[refs][di][dik]['pos'][pi]['ad_ct'] == 'na':
                                        if "AD" in condict.keys():
                                            full_dict[refs][di][dik]['pos'][pi]['ad_ct'] = condict["AD"]
                                    if full_dict[refs][di][dik]['pos'][pi]['con_ct'] == 'na':
                                        if "Control" in condict.keys():
                                            full_dict[refs][di][dik]['pos'][pi]['con_ct'] = condict["Control"]
                                    
                                    query2 = f"SELECT DISTINCT edit_dist FROM {table} WHERE ref_base = ? AND a1type = ? and poss = ?"
                                    que2 = co.execute(query2, (refs,di,pos)).fetchall()
                                    global_dist = que2[0][0]
                                    if not 'global_dist' in full_dict[refs][di][dik]['pos'][pi].keys():
                                        full_dict[refs][di][dik]['pos'][pi]['global_dist'] = global_dist
                                    
                                    query3 = f"SELECT DISTINCT con,rids,a1prop FROM {table} WHERE ref_base = ? AND poss = ?"
                                    que3 = co.execute(query3, (refs,pos)).fetchall()
                                    for q3 in que3:
                                        condis,props = q3[0],q3[2]
                                        if condis == "AD":
                                            ad_list.append(props)

                                        if condis == "Control":
                                            con_list.append(props)
                                    
                
                    prop_diff_pval = 'na'
                    higher_prop_con = 'na'
                    if con_list and ad_list:
                        stat, p_value = mannwhitneyu(ad_list, con_list, alternative='two-sided')
                        prop_diff_pval = p_value
                        ad_sum  = sum(ad_list)
                        con_sum = sum(con_list)
                        ad_prop,con_prop = ad_sum/len(ad_list),con_sum/len(con_list)
                        if ad_prop > con_prop:
                            higher_prop_con = "AD"
                        if con_prop > ad_prop:
                            higher_prop_con = "Control"
                    full_dict[refs][di][dik]['dif_val'] = prop_diff_pval
                    full_dict[refs][di][dik]['dif_higher'] = higher_prop_con
                    
                    for re in full_dict.keys():
                        for typing in full_dict[re].keys():
                            for loc in full_dict[re][typing].keys():
                                diff_pval,diff_high = full_dict[re][typing][loc]['dif_val'],full_dict[re][typing][loc]['dif_higher']
                                for position in full_dict[re][typing][loc]['pos'].keys():
                                    ads,cons,glob,hgv = full_dict[re][typing][loc]['pos'][position]['ad_ct'],full_dict[re][typing][loc]['pos'][position]['con_ct'],full_dict[re][typing][loc]['pos'][position]['global_dist'],full_dict[re][typing][loc]['pos'][position]['hgvs']
                                    outrow = [re,position,ads,cons,typing,glob,loc,hgv,diff_pval,diff_high]
                                    final_write.writerow(outrow)
        print(refs)
    return

def get_sample_dists(dbp):
    

    outfile = open('global_sample_dists.csv','w')
    outwrite = csv.writer(outfile)
    header = ['global_dist','var_type','ad_samps','con_samps']
    outwrite.writerow(header)
    con = sql.connect(db)
    co  = con.cursor()
    tab = "All_Info_TopGenes_tab"
    info_dict = {}
    samps = co.execute(f"SELECT DISTINCT edit_dist  COUNT(DISTINCT rids) FROM {tab} GROUP BY ").fetchall()
    return

if __name__=='__main__':
    dbp = sys.argv[1]
    ad_only = sys.argv[2]
    ser = sys.argv[3]
    nov = sys.argv[4]
    filt_loc = sys.argv[5]
    #dist_cts = get_edit_dist_cts(dbp)
    #samples = get_sample_info(dbp)
    #edit_types = check_distinct_by_edit_type(dbp)
    #complexes = table2_info(dbp)
    #global_local = global_local_edit_types(dbp)
    #filts_deets = local_set_deets(dbp,filt_loc)
    #per_gene_investigate(dbp,filt_loc)
    pcas = create_dat_for_PCA(dbp)
