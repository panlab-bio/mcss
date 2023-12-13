# -*- coding: utf-8 -*-
import pandas as pd
import sys
import os
import pickle
path_sp = sys.argv[1]
path_sim = sys.argv[2]#我直接输入path_sim原版就好了啊
path_multi = int(sys.argv[3])
# print(path_multi)
if path_multi==0:
    # print("一个输入")
    df_sp = pd.read_csv(path_sp,header=None,sep="\t")
    df_sp.columns=["Species"]
    name_sp_list = [df_sp["Species"].to_list()]
    print(name_sp_list)
else:
    # print("多组输入")
    sp_list = sorted(os.listdir(path_sp))
    name_sp_list = []
    for sp_a in sp_list:
        path_sp_a_sample = os.path.join(path_sp,sp_a)
        # print(path_sp_a_sample)
        df_sp = pd.read_csv(path_sp_a_sample,header=None,sep="\t")
        df_sp.columns=["Species"]
        name_sp_list.append(list(df_sp["Species"].to_list()))

# if path_sim.endswith("/"):
#     path_sim_root = path_sim.replace("community/","")
# else:
#     path_sim_root = path_sim.replace("community","")
path_save = os.path.join(path_sim,"species","sp_list.pkl")

with open(path_save,"wb") as f:
    pickle.dump(name_sp_list,f)
df_name_sp_list = pd.DataFrame(name_sp_list)
df_name_sp_list.columns=["species"]
path_save_df = path_save.replace("pkl","csv")
df_name_sp_list.to_csv(path_save_df,index=False)
