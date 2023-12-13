# -*- coding: utf-8 -*-
#生成模拟的树
#输入的是sim路径 脚本路径 物种个数 env
import numpy as np
import pandas as pd
import itertools
import subprocess
import sys
import os
import random
import pickle 
import copy
from lib.dis_generate import *

path_sim = sys.argv[1]#这里的path_sim也不是community
path_abs = sys.argv[2]
cnt_sample = int(sys.argv[3])
env = sys.argv[4]
trim_len = sys.argv[5]
path_tax = sys.argv[6]
flag_acc_int = int(sys.argv[7])

if flag_acc_int == 0:
    flag_dis = True
else:
    flag_dis = False

path_tree =  os.path.join(path_abs,"data/n_ary_tree")
path_tree_name =  os.path.join(path_tree,"data_tree_name.pkl")
path_tree_dict =  os.path.join(path_tree,"dis_all_dict_alltree.pkl")
path_tree_data =  os.path.join(path_tree,"data_sp_alltree.pkl")
path_tree_dis =  os.path.join(path_tree,"dis_tree.pkl")
path_tree_shift =  os.path.join(path_tree,"shift_tree.pkl")
# print(path_tree_shift)

with open(path_tree_data,"rb") as f:
    data_alltree = pickle.load(f)  
# 总树的字典
with open(path_tree_dict,"rb") as f:
    dict_tree = pickle.load(f)
    
with open(path_tree_name,"rb") as f:
    datad_name = pickle.load(f)
    
with open(path_tree_dis,"rb") as f:
    dis_tree = pickle.load(f)
with open(path_tree_shift,"rb") as f:
    shift_tree = pickle.load(f)

path_train = os.path.join(path_abs,"data/env",env)
path_name = os.path.join(path_train,"name_env.pkl")
path_data = os.path.join(path_train,"data_env.pkl")
path_dis = os.path.join(path_train,"dis_env.pkl")
path_shift = os.path.join(path_train,"shift_env.pkl")
# path_dict = os.path.join(path_train,"dict_env.pkl")


# with open(path_dict,"rb") as f:
#      dict_tree_train = pickle.load(f)  

with open(path_data,"rb") as f:
    data_tree_train = pickle.load(f)
    
with open(path_name,"rb") as f:
    name_tree_train = pickle.load(f)
    
with open(path_dis,"rb") as f:
    dis_tree_train = pickle.load(f)
    
with open(path_shift,"rb") as f:
    shift_tree_train = pickle.load(f)
    
taxonomy = pd.read_csv(path_tax,sep="\t")
list_ar = taxonomy[taxonomy["Domain"]=="Archaea"]["Phylum"].to_list()
list_bac = taxonomy[taxonomy["Domain"]=="Bacteria"]["Phylum"].to_list()

name_phy_bac = []
name_phy_ar = []
position_phy = dict()
for i,name_tree in enumerate(name_tree_train):
    for name_p_one in name_tree[1]:
        if name_p_one in list_ar:
            name_phy_ar.append(name_p_one)
        else:
            name_phy_bac.append(name_p_one)
    for sp_i,name_sp in enumerate(name_tree[1]):
        if name_sp not in  position_phy.keys():
            position_phy[name_sp]= [(i,sp_i,data_tree_train[i][1][sp_i])]
        else:
            position_phy[name_sp].append((i,sp_i,data_tree_train[i][1][sp_i]))
        
uni_name_phy_bac,cnt_name_phy_bac = np.unique(name_phy_bac,return_counts=True)
cnt_name_phy_std_bac = [ cnp/sum(cnt_name_phy_bac) for cnp in cnt_name_phy_bac]
# print("bac",len(name_phy_bac),len(uni_name_phy_bac),sum(cnt_name_phy_std_bac))

uni_name_phy_ar,cnt_name_phy_ar = np.unique(name_phy_ar,return_counts=True)
cnt_name_phy_std_ar = [ cnp/sum(cnt_name_phy_ar) for cnp in cnt_name_phy_ar]
# print("ar",len(name_phy_ar),len(cnt_name_phy_std_bac),sum(cnt_name_phy_std_bac))

#最大最小的物种个数
cnt_sp_70 = []

for i in range(len(data_tree_train)):
    cnt_sp_70.append(len(name_tree_train[i][-1]))    
mean_cnt = np.mean(cnt_sp_70)
std_cnt = np.std(cnt_sp_70)
min_cnt = mean_cnt - 3*std_cnt
max_cnt = mean_cnt + 3*std_cnt

if trim_len==0:#
    min_cnt = 0
    max_cnt = 1000

path_meanshift = os.path.join(path_abs,"data/meanshift",env)

data_tree_root,dis_tree_root,dis_shift_tree_root = get_root_node(data_alltree,dis_tree,shift_tree) 
datad_name_root = [["root"]]+datad_name
name_sp_list = []
# len_sp_list = []
# list_dis2domian_sim =[]
# list_dl1 = []
# list_dl2 = []
cnt_sub_tree = 0 #记录得到的子树的个数

while cnt_sub_tree<cnt_sample:
    # cnt_sub_tree+=1
    data_sim,dis_sp_sim,dis_shift_sim = get_sp_cnt_data_level(path_meanshift,
    data_tree_train,dis_tree_train,shift_tree_train,position_phy,
    uni_name_phy_bac,cnt_name_phy_std_bac,uni_name_phy_ar,cnt_name_phy_std_ar)
    
    if sum(data_sim[-1])<=max_cnt and sum(data_sim[-1])>=min_cnt:
        # print("---",min_cnt,max_cnt,sum(data_sim[-1]))
        cnt_sub_tree+=1
        
        data_sim_root,dis_sim_root,dis_shift_sim_root = get_root_node(data_sim,dis_sp_sim,dis_shift_sim)  
        flag,_ ,_= is_subtree2(data_sim_root,data_tree_root)

        if flag:
            dis_12,name_sp_12 = get_tree_dis_new_3(data_sim_root,dis_sim_root,
                    data_tree_root,dis_tree_root,datad_name_root,0,False,flag_dis)
            name_sp_list.append(list(name_sp_12[-1]))
        print("generate ",len(name_sp_12[-1]),"species","dis",dis_12)
            
            


path_save = os.path.join(path_sim,"species","sp_list.pkl")

with open(path_save,"wb") as f:
    pickle.dump(name_sp_list,f)
df_name_sp_list = pd.DataFrame(name_sp_list)
df_name_sp_list.columns=["species"]
path_save_df = path_save.replace("pkl","csv")
df_name_sp_list.to_csv(path_save_df,index=False)
