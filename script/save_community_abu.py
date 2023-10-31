# -*- coding: utf-8 -*-
#将community模式下的输入改成pkl文件
#这里用户自己输入的可能不是pkl
import numpy as np
import pandas as pd
import sys
import os
import pickle 

path_sp_root = sys.argv[1]
path_abu_root = sys.argv[2]
multi = int(sys.argv[3])
path_sim = sys.argv[4]
path_save = os.path.join(path_sim,"community/abu/abu_user.pkl")
if multi==0:
    #一个文件
    abu_sample_sim_all = []
    path_abu = path_abu_root
    abu_real = pd.read_csv(path_abu,header=None,sep="\t")
    abu_real.columns=["abu"]
    abu_sample_sim = abu_real["abu"].to_list()
    abu_sample_sim_all.append(list(abu_sample_sim))
    
    with open(path_save,"wb") as f:
        pickle.dump(abu_sample_sim_all,f)
else:
    list_sp = sorted(os.listdir(path_sp_root))
    abu_sample_sim_all = []
    for sp in list_sp:
        name_abu = sp.replace("sp","abu")
        path_abu = os.path.join(path_abu_root,name_abu)
        abu_real = pd.read_csv(path_abu,header=None,sep="\t")
        abu_real.columns=["abu"]
        abu_sample_sim = abu_real["abu"].to_list()
        abu_sample_sim_all.append(list(abu_sample_sim))
    with open(path_save,"wb") as f:
        pickle.dump(abu_sample_sim_all,f)

