# -*- coding: utf-8 -*-
# Generate the simulated tree based on community.
import pandas as pd
import sys
import os
# import pickle
import json
path_sp = sys.argv[1]
path_sim = sys.argv[2]
path_multi = int(sys.argv[3])

if path_multi==0: #One input.
    
    df_sp = pd.read_csv(path_sp,header=None,sep="\t")
    df_sp.columns=["Species"]
    name_sp_list = [df_sp["Species"].to_list()]
    print(name_sp_list)
else: # Multiple inputs.
    
    sp_list = sorted(os.listdir(path_sp))
    name_sp_list = []
    for sp_a in sp_list:
        path_sp_a_sample = os.path.join(path_sp,sp_a)
        
        df_sp = pd.read_csv(path_sp_a_sample,header=None,sep="\t")
        df_sp.columns=["Species"]
        name_sp_list.append(list(df_sp["Species"].to_list()))


path_save = os.path.join(path_sim,"species","sp_list.json")

with open(path_save,"w") as f:
    json.dump(name_sp_list,f)
for nsl_i,nsl in enumerate(name_sp_list):
    path_save_df = path_save.replace(".json","_"+str(nsl_i)+".csv") 

    df_name_sp_list = pd.DataFrame(nsl)
    df_name_sp_list.columns=["species"]
    df_name_sp_list.to_csv(path_save_df,index=False)
