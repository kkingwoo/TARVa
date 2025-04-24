import csv
import pandas as pd

class ParseVep:
    @staticmethod
    def parse_v(vfile,set_name):
        g_dict = {}
        v = pd.read_csv(vfile,delimiter='\t',usecols=['#Uploaded_variation','SYMBOL','Gene','PolyPhen','Location'],low_memory=False)
        cols = ['#Uploaded_variation','SYMBOL','Gene','PolyPhen','Location']
        infos = v[cols].values.tolist()
        for rows in infos:
            g,hu = rows[2],rows[1]
            if not g in g_dict.keys():
                g_dict[g] = hu
        return set_name,g_dict
