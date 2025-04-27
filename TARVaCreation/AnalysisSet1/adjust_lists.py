class ListLens:
    @staticmethod
    def adjust_lens(listing,ct,out_dict,key):
        if len(listing) < ct:
            zeros = ct - len(listing)
            listing += [0] * zeros    
        out_dict[key] = []
        out_dict[key].append(listing)
        
        return out_dict

    def alg_bytype(df,con,type_key,type_value,percent_key):
        new_df = df[(df['condition'] == con) & (df[type_key] == type_value)][['rid', percent_key]].values.tolist()
        
        return new_df

    def eag_by_position(ads,controls,gene):
        single_calls,single_condition_mult_calls,both_conditions = [],[],[]
        if not gene in controls.keys():
            for position in ads[gene].keys():
                count = ads[gene][position]
                if count == 1:
                    single_calls.append(position)
                else:
                    single_condition_mult_calls.append(position)
        else:
            for position in ads[gene].keys():
                ad_ct = ads[gene][position]
                if not position in controls[gene].keys():
                    if ad_ct == 1:
                        single_calls.append(position)
                    else:
                        single_condition_mult_calls.append(position)
                else:
                    both_conditions.append(position)
            for position in controls[gene].keys():
                con_ct = controls[gene][position]
                if not position in ads[gene].keys():
                    if con_ct == 1:
                        single_calls.append(position)
                    else:
                        single_condition_mult_calls.append(position)
                else:
                    both_conditions.append(position)

        return gene,single_calls,single_condition_mult_calls,both_conditions

    def get_data(long_list,append_list):
        for l in long_list:
            append_list.append(l[1])
        return append_list

