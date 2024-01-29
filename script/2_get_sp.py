# -*- coding: utf-8 -*-
# Process the  .out and .rep files
# Truncate species based on their abundance.


import numpy as np
import pandas as pd
import os
import sys

path_data = sys.argv[1]
path_s1 = sys.argv[2]
abu_sp = float(sys.argv[3])
start = int(sys.argv[4])
end = int(sys.argv[5])


path_abundance_root = os.path.join(path_data,"kraken/abundance/")
path_read_cnt =  os.path.join(path_data,"kraken/read_cnt/read_cnt.txt")
path_cut = os.path.join(path_data,"sp_cut")

# get species of the sample 
def get_sp_list(path_read_cnt,path_abundance_root,path_cut,path_s1 = path_s1):
    
    
    s1_sp = pd.read_csv(path_s1,header=None,sep="\t",dtype=str)
    s1_sp.columns = ["name","s1"] 
    # Obtain the correspondence between sp and s1, then merge with each sample, and save the merged results in the sp_cut folder.
    abundance_name = sorted(os.listdir(path_abundance_root))[start:end]
    
    for i in abundance_name:
        name = i.replace("_abundance.txt","_cutsp.txt")
        path_name = os.path.join(path_cut,name)
        path_abu = os.path.join(path_abundance_root,i)
        
        rep = pd.read_csv(path_abu,sep="\t",header=None)
        rep.columns =["name","abu"]
        list_sp_cnt = rep["abu"].to_list()
        read_cnt_flag = int(sum(list_sp_cnt)*abu_sp) 
        
        if read_cnt_flag>0:

            mask = (rep["abu"]<=read_cnt_flag) 
            rep_drop = rep.drop(rep[mask].index)
            df_merge  = pd.merge(rep_drop,s1_sp,on=["name"])
            df_merge_c2 = df_merge.loc[:,["name","s1"]]
            print("sp",len(df_merge_c2))
            len_merge = len(df_merge_c2)
            
            df_merge_c2.to_csv(path_name,sep="\t",index=None,header=None)

get_sp_list(path_read_cnt,path_abundance_root,path_cut) 
