import sqlite3 as sql
import pandas as pd
import numpy as np
import csv
from scipy.stats import mannwhitneyu

class EditTypes:
    @staticmethod
    def edits(d,g):
        transcript_ct = 0
        type_list = ['A-I','A-C','A-T','C-U','C-G','C-A','G-A','G-T','G-C','T-A','T-C','T-G','ins','del']
        db_con = sql.connect(d,check_same_thread=False)
        db = db_con.cursor()
        check = []
        trans_dict = {}
        s_tab = 'sample_tab'
        g_tab = 'gtf_tab'
        transcripts = f"SELECT DISTINCT ref_id,tid FROM {g_tab} where ref_id = ? GROUP BY ref_id,tid"
        trans = db.execute(transcripts,(g,)).fetchall()
        full_pos_list = []
        for tra in trans:
            lens_query = f"SELECT DISTINCT tid,trans_start,trans_end FROM {g_tab} WHERE tid= ? GROUP BY tid"
            lens = db.execute(lens_query, (tra[1],)).fetchall()
            for l in lens:
                transcript_ct+=1
                trr = l[0]
                trans_dict[trr] = ''
                start,end = l[1],l[2]
                trans_dict[trr] = (start,end)
                full_pos_list.append(start)
                full_pos_list.append(end)
        full_pos_list.sort()
        full_set = list(set(full_pos_list))
        full_set.sort()

        g_start,g_end = full_set[0],full_set[-1]
        gene_len = g_end-g_start
        gene_pos_wt = 1/gene_len
        
        hgvs_list = []

        dist_dict = {}

        unique_positions_dict = {}
        unique_positions_dict['+'] = {}
        unique_positions_dict['-'] = {}
        
        props_dict  = {}
        props_dict['AD'] = {}
        props_dict['Control'] = {}
        for tt in type_list:
            props_dict['AD'][tt] = {}
            props_dict['Control'][tt] = {}
        
        main_out = []

        neg_strand_call,neg_strand_rna = ('A','C','T','G'),('T','G','A','C')
        
        from_samps = f"SELECT DISTINCT condition,rid,pos,ref_id,a1_ad,a2_ad,dp,alt1,alt2,ref,strand,trans_pos,chrom,tid FROM {s_tab} WHERE ref_id = ? GROUP BY condition,rid,pos,ref_id,a1_ad,a2_ad,dp"      
        samps = db.execute(from_samps,(g,)).fetchall()
        for s in samps:
        
            con,rids,poss,refid,ref_base,a1ad,a2ad,dp,alt1,alt2,ref,strand,trans_pos,chrom,tid= s[0],s[1],s[2],s[3],s[3].split('_')[0],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13]
            
            pos = str(poss)


            ##distinct positions and information for checking against REDi DB
            if not chrom in unique_positions_dict[strand].keys():
                unique_positions_dict[strand][chrom] = []
            if not pos in unique_positions_dict[strand][chrom]:
                unique_positions_dict[strand][chrom].append(pos)
            
            ##create a dictionary for determining which distribution type the gene belongs in  Does not double count for same positions onopposite strands...
            if not pos in dist_dict.keys():
                dist_dict[pos] = {}
            if not con in dist_dict[pos].keys():
                dist_dict[pos][con] = []
            if not rids in dist_dict[pos][con]:
                dist_dict[pos][con].append(rids)

            
            hgvs = ''
            a1type = ''
            
            

            ## Assign edit-type based on strand, ref, alt1 information
            if strand == '+':
                if alt1 == '*':
                    a1type = 'del'
                elif len(ref) > len(alt1):
                    a1type = 'del'
                elif len(ref) < len(alt1):
                    a1type = 'ins'
                elif ref == 'A' and alt1 == 'G':
                    a1type = 'A-I'
                elif ref == 'C' and alt1 == 'T':
                    a1type = 'C-U'
                else:
                    a1type = f"{ref}-{alt1}" if len(ref) == len(alt1) else 'unknown'
            elif strand == '-':
                ri = neg_strand_call.index(ref) if ref in neg_strand_call else -1
                ai = neg_strand_call.index(alt1) if alt1 != '*' and alt1 in neg_strand_call else -1
                ro,ao = (neg_strand_rna[ri], neg_strand_rna[ai]) if ri != -1 and ai != -1 else ('N', 'N')
                if alt1 == '*':
                    a1type = 'del'
                elif len(ref) > len(alt1):
                    a1type = 'del'
                elif len(ref) < len(alt1):
                    a1type = 'ins'
                elif ref == 'T' and alt1 == 'C':
                    a1type = 'A-I'
                elif ref == 'G' and alt1 == 'A':    
                    a1type = 'C-U'
                else:
                    a1type = f"{ro}-{ao}" if len(ref) == len(alt1) else 'unkown'        

           ## Remove hgvs throughout script... not relevant to study 
            ty = type_list.index(a1type)
            enst = tid.replace('_','.')
            if ty in range(0,14):
                hgvs = f"{enst}:c.{trans_pos}{ref}>{alt1}"
            elif ty == 14:
                po_in_end = trans_pos+len(alt1)
                hgvs = f"{enst}:c.{trans_pos}_{po_in_end}ins{alt1}"
            else:
                po_in_end = trans_pos+len(ref)
                hgvs = f"{enst}:c.{trans_pos}_{po_in_end}del"
            if hgvs not in hgvs_list:
                hgvs_list.append(hgvs)
            
            ##Create dictionary for running mwu by edit-type between the two conditions, within the same gene
            if rids not in props_dict[con][a1type].keys():
                props_dict[con][a1type][rids]=float(0.0)
            a1prop = float((a1ad/dp)*gene_pos_wt)
            props_dict[con][a1type][rids]+=a1prop   

            check.append((con,rids,poss,refid,ref_base,a1ad,dp,ref,alt1,a1type,strand,chrom,hgvs,a1prop))

        props_cts = []
        for key in props_dict.keys():
            tup = tuple()
            for k in props_dict[key].keys():
                props_list = []
                for sampling in props_dict[key][k].keys():
                    amt = props_dict[key][k][sampling]
                    props_list.append(amt)
                tup = (key,k,props_list)
                props_cts.append(tup)
        ad_list,cont_list = [],[]
        for tu in props_cts:
            con,types,lists = tu[0],tu[1],tu[2]
            if con == "AD":
                add_zeros = 29 - len(lists)
                lists.extend([float(0.0)]*add_zeros)
                ad_list.append((types,lists))
            if con == "Control":
                add_zeros = 30 - len(lists)
                lists.extend([float(0.0)]*add_zeros)
                cont_list.append((types,lists))

        out_stats = []
        for infos in range(0,len(ad_list)):
            typing = ad_list[infos][0]
            ad,conn = ad_list[infos][1],cont_list[infos][1]
            stat,p_value = mannwhitneyu(ad,conn)
            if float(p_value) <= float(0.05):
                out_stats.append((refid,typing,stat,p_value,hgvs,ad,conn))


        check_edit_dist = []
        check_outs = {}
        check_outs['AD'] = {}
        check_outs['Control'] = {}

        disted = ['unique','common']

        for i in range(0,len(disted)):
            di = disted[i]
            check_outs['AD'][di] = 0
            check_outs['Control'][di] = 0
        for positions in dist_dict.keys():
            for conditions in dist_dict[positions].keys():
                if dist_dict[positions][conditions]:
                    sampsize = len(dist_dict[positions][conditions])
                    check_tupe = (positions,conditions,sampsize)
                    check_edit_dist.append(check_tupe)
        for checking in check_edit_dist:
            co,num = checking[1],checking[2]
            if num == 1:
                check_outs[co]['unique']+=num
            else:
                check_outs[co]['common']+=num

        gene_dist_type = ''

        ad_unique,ad_common = check_outs['AD']['unique'],check_outs['AD']['common']
        control_unique,control_common = check_outs['Control']['unique'],check_outs['Control']['common']

        unique_sums,common_sums = ad_unique+control_unique,ad_common+control_common

        ad_common_perc,control_common_perc,ad_unique_perc,control_unique_perc = 0,0,0,0

        if common_sums != 0:
            ad_common_perc = float(ad_common/common_sums)
            control_common_perc = float(control_common/common_sums)
        if unique_sums != 0:
            ad_unique_perc = float(ad_unique/unique_sums)
            control_unique_perc = float(control_unique/unique_sums)

        if (ad_unique+ad_common) >= 1 and (control_unique+control_common) == 0:
            gene_dist_type = 'AD_only'
        elif (control_unique+control_common) >= 1 and (ad_unique+ad_common) == 0:
            gene_dist_type = 'Control_only'
        elif (ad_unique_perc+ad_common_perc) >= float(0.75):
            gene_dist_type = 'Primarily_AD'
        elif (control_unique_perc+control_common_perc) >= float(0.75):
            gene_dist_type = 'Primarily_Control'
        else:
            gene_dist_type = 'Both'

        out_string = f"Gene {ref_base} has {transcript_ct} possible transcripts associated\n"
        return refid,gene_dist_type,out_stats,unique_positions_dict,check,hgvs_list,out_string,dist_dict
