import os
import pandas as pd
import numpy as np

class ReadMstrg:
    @staticmethod
    def read_string_tpm(con,stri,samp,od):
        out_list = []
        path = stri+samp+'/t_data.ctab'
        tab = pd.read_csv(path,delimiter='\t')
        total_fpkm = tab['FPKM'].sum()
        for refs in od.keys():
            for tid in od[refs]:
                ti = tid.replace('_','.')
                infos = tab.loc[tab['t_name']== ti,['length','FPKM']].values.tolist()
                for i in infos:
                    tpm = float(0.0)
                    length,fpkm = i[0],i[1]
                    if total_fpkm > float(0.0):
                        tpm = float((fpkm / total_fpkm) * 1e6)
                    tupe = (con,refs,ti,length,tpm)
                    out_list.append(tupe)
        return out_list

    @staticmethod
    def read_string_fpkm(con,stri,samp,od):
        out_list = []

        path = stri+samp+'/t_data.ctab'
        tab = pd.read_csv(path,delimiter='\t')
        for refs in od.keys():
            for tid in od[refs]:
                outs = [refs,tid,samp,con]
                ti = tid.replace('_','.')
                infos = tab.loc[tab['t_name']== ti,'FPKM'].values.tolist()
                for i in infos:
                    if not float(i) == float(0.00):
                        outs.append(float(i))
                        out_list.append(outs)

        return out_list
