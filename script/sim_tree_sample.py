# -*- coding: utf-8 -*-
#把三十个测试数据的data_test,dis_test,shift_test都保存起来
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
from scipy.stats import lognorm
# from distfit import distfit

path_sim = sys.argv[1]
path_tax = sys.argv[2]
path_nw = sys.argv[3]
path_abs = sys.argv[4]
flag_acc_int = int(sys.argv[5])
cnt_sample =  int(sys.argv[6])
same_env = int(sys.argv[7])


if flag_acc_int == 0:
    flag_dis = True
else:
    flag_dis = False

# if change_len_sp == 0:#0的时候不change
#     flag_len_sp = False
#     print("不随机确定")
# else:
#     flag_len_sp = True
    
if same_env == 0:#0的时候不change
    flag_env = False
    # print("不是同一环境")
else:
    flag_env = True
    
path_tree =  os.path.join(path_abs,"data/n_ary_tree")
path_tree_dict =  os.path.join(path_tree,"dis_all_dict_alltree.pkl")
with open(path_tree_dict,"rb") as f:
    dict_tree = pickle.load(f)
    
data_tree_30 =[]
dis_tree_30 =[]
shift_tree_30 =[]
name_tree_30 = []
taxonomy=pd.read_csv(path_tax,sep="\t")
path_name_root = os.path.join(path_sim,"data/sp_cut/")
path_com = os.path.join(path_sim,"data/sp_cut/")
list_name = sorted(os.listdir(path_name_root))

def get_min_sp_list_dis(dis_sp_domain,dis_list):
    min_dis = 100
    
    for dis in dis_list:
        if abs(dis_sp_domain-dis)<min_dis:
            min_dis = abs(dis_sp_domain-dis)
    return min_dis

if len(list_name)==1:#用户只提供了单个样本
    flag_one = True
else:
    flag_one = False
list_sp_genus_ori_all = [] #记录每个样本的物种池
list_abu_all = []
path_abu_root = os.path.join(path_sim,"data/kraken/abundance")
for name in list_name:
    # print(name)\
    abu = name.replace("cutsp","abundance")
    path_abu = os.path.join(path_abu_root,abu)
    df_abu = pd.read_csv(path_abu,header=None,sep="\t")
    df_abu.columns=["Species","cnt"]
    
    path_name = os.path.join(path_name_root,name)
    # print(path_name)
    sp_cut=pd.read_csv(path_name,header=None,usecols=[0],sep="\t")
    sp_cut.columns=["Species"]
    merge_name = pd.merge(taxonomy,sp_cut,on=["Species"])
    df_merge = pd.merge(sp_cut,df_abu,on=["Species"])
    cnt_sample_list = df_merge["cnt"].to_list()
    cnt_sum = sum(cnt_sample_list)
    abu_sample = [float("{:.4f}".format(cs/cnt_sum)) for cs in cnt_sample_list]
    if flag_env:
        print("abu env")
        list_abu_all+=abu_sample #一个环境多个样本的时候直接把所有的abu加起来
    else:#如果一个环境一个样本
        print("abu distribution")
        # path_abu_save = os.path.join(path_sim,"abu/abu.pkl")
        # print(len(abu_sample),abu_sample)
        # dist = distfit(distr='lognorm',verbose = 0)
        # x_abu_sample = np.array(abu_sample)
        # dist.fit_transform(x_abu_sample)
        # list_abu_all.append(dist)
        x_abu_sample = np.array(abu_sample)
        shape_ln, loc_ln, scale_ln = lognorm.fit(x_abu_sample, floc=0)
        dict_lognorm = {"shape":shape_ln,"loc":loc_ln,"scale":scale_ln}
        list_abu_all.append(dict_lognorm)
        # print("abu sample")
        
        
    dis_sp_list = []#得到了原始的距离列表
    for sp_genus in merge_name["Species"].to_list():
        if sp_genus in dict_tree["d"][0]["Archaea"][-1].keys():
            dis_sp_domain = dict_tree["d"][0]["Archaea"][-1][sp_genus]
        else:
            dis_sp_domain = dict_tree["d"][0]["Bacteria"][-1][sp_genus]
        dis_sp_list+=[dis_sp_domain]
        
#     if flag_len_sp :#长度可变的情况下，我们需要随机选择物种,如果是单环境多样本，不用考虑这个
#         # len_sp_sample = len(merge_name)
#         # len_sp_new_list = list(range(round(len_sp_sample-0.1*len_sp_sample),round(len_sp_sample+0.1*len_sp_sample)))
#         # len_sp_new = random.choice(len_sp_new_list)
#         merge_name_genus = merge_name["Genus"].unique()
#         taxonomy_genus = taxonomy[taxonomy["Genus"].isin(merge_name_genus)].reset_index(drop=True)
        
#         list_sp_genus_ori =[]
#         for sp_genus in taxonomy_genus["Species"].to_list():
#             if sp_genus in dict_tree["d"][0]["Archaea"][-1].keys():
#                 dis_sp_domain = dict_tree["d"][0]["Archaea"][-1][sp_genus]
#             else:
#                 dis_sp_domain = dict_tree["d"][0]["Bacteria"][-1][sp_genus]
#             min_dis_sp = get_min_sp_list_dis(dis_sp_domain,dis_sp_list)
#             if min_dis_sp<=0.1:
#                 list_sp_genus_ori.append(sp_genus)
        
#         list_sp_genus_ori_all.append(list_sp_genus_ori)        
        # list_sp_genus = random.sample(list_sp_genus_ori,len_sp_new)
        # merge_name_sp = taxonomy[taxonomy["Species"].isin(list_sp_genus)].reset_index(drop=True)
    # else:
    merge_name_sp = copy.deepcopy(merge_name)
        
    data_sp_real,datad_name30 = get_num_new_2(merge_name_sp)
    dis_all_dict_real,data_sp_real_2 = get_dis_list_alltree(merge_name_sp,path_tax,path_nw,path_sim,path_abs) 
    #这个得到的data_sp_real好像不对,为了省事让他是空
    dis_tree30,shift_tree30 = get_dis_tree(dis_all_dict_real,datad_name30,data_sp_real,merge_name_sp)
    data_tree_30.append(data_sp_real)
    dis_tree_30.append(dis_tree30)
    shift_tree_30.append(shift_tree30)
    name_tree_30.append(datad_name30)
if flag_env:
    abu_name = "abu_env.pkl"
else:
    abu_name = "abu_sample.pkl"
    # path_abu_save = os.path.join(path_sim,"abu/abu.pkl")#这个时候abu.pkl存的不是分布，而是物种丰度列表的和
    # with open(path_abu_save) as f:
    #     pickle.dump(list_abu_all,f)
path_abu_save = os.path.join(path_sim,"abu",abu_name)#这个时候abu.pkl存的不是分布，而是物种丰度列表的和
print(path_abu_save)
with open(path_abu_save,"wb") as f:
    pickle.dump(list_abu_all,f)
# with open(path_abu_save,"rb") as f:
#     list_abu_all = pickle.load(f)
# dist = list_abu_all[0]
# print("dist",dist.generate(1))
# else:
#     path_abu_save = os.path.join(path_sim,"abu/abu.pkl")
#     dist.save(path_abu_save)
#不用保存了，因为在这个文件里自己调用了
#         path_sample = os.path.join(path_sim,"data/sample")
#         path_name = os.path.join(path_sample,"name_sample.pkl")
#         path_data = os.path.join(path_sample,"data_sample.pkl")
#         path_dis = os.path.join(path_sample,"dis_sample.pkl")
#         path_shift = os.path.join(path_sample,"shift_sample.pkl")
#         path_dict = os.path.join(path_sample,"dict_sample.pkl")

#         with open(path_dict,"wb") as f:
#             pickle.dump(dis_all_dict_real,f)

#         with open(path_name,"wb") as f:
#             pickle.dump(name_tree_30,f)

#         with open(path_data,"wb") as f:
#             pickle.dump(data_tree_30,f)

#         with open(path_dis,"wb") as f:
#             pickle.dump(dis_tree_30,f)

#         with open(path_shift,"wb") as f:
#             pickle.dump(shift_tree_30,f)


path_tree_name =  os.path.join(path_tree,"data_tree_name.pkl")
path_tree_data =  os.path.join(path_tree,"data_sp_alltree.pkl")
path_tree_dis =  os.path.join(path_tree,"dis_tree.pkl")
path_tree_shift =  os.path.join(path_tree,"shift_tree.pkl")
# print(path_tree_shift)

with open(path_tree_data,"rb") as f:
    data_alltree = pickle.load(f)  
    
with open(path_tree_name,"rb") as f:
    datad_name = pickle.load(f)
    
with open(path_tree_dis,"rb") as f:
    dis_tree = pickle.load(f)
with open(path_tree_shift,"rb") as f:
    shift_tree = pickle.load(f)

data_tree_root,dis_tree_root,dis_shift_tree_root = get_root_node(data_alltree,dis_tree,shift_tree)  
datad_name_root = [["root"]]+datad_name
name_sp_list = []
if flag_env:#就是一个环境多个样本
    #这和env模式类似，但是我们不用mean_shift了，直接抽样
    list_ar = taxonomy[taxonomy["Domain"]=="Archaea"]["Phylum"].to_list()
    list_bac = taxonomy[taxonomy["Domain"]=="Bacteria"]["Phylum"].to_list()
    name_phy_bac = []
    name_phy_ar = []
    position_phy = dict()
    for i,name_tree in enumerate(name_tree_30):
        for sp_i,name_sp in enumerate(name_tree[1]):
            if name_sp in list_ar:
                name_phy_ar.append(name_sp)
            elif name_sp in list_bac:
                name_phy_bac.append(name_sp)
            else:
                print("error,not ar or bac")
            
            if name_sp not in  position_phy.keys():
                position_phy[name_sp]= [(i,sp_i,name_tree_30[i][1][sp_i])]
            else:
                position_phy[name_sp].append((i,sp_i,name_tree_30[i][1][sp_i]))
                
    uni_name_phy_bac,cnt_name_phy_bac = np.unique(name_phy_bac,return_counts=True)
    cnt_name_phy_std_bac = [ cnp/sum(cnt_name_phy_bac) for cnp in cnt_name_phy_bac]
    uni_name_phy_ar,cnt_name_phy_ar = np.unique(name_phy_ar,return_counts=True)
    cnt_name_phy_std_ar = [ cnp/sum(cnt_name_phy_ar) for cnp in cnt_name_phy_ar]
    # print("uni_name_phy_bac,uni_name_phy_ar",uni_name_phy_bac,uni_name_phy_ar)
    len_d_list =[[],[]]#domian count
    # len_p_list = []
    name_p_list = []                
    for dt30,nt30 in zip(data_tree_30,name_tree_30):
        df_name_sp = pd.DataFrame(nt30[1],columns=["Phylum"])
        # print(df_name_sp)
        df_merge_p = pd.merge(df_name_sp,taxonomy,on=["Phylum"])
        #0 archaea 1 bacteria
        l_archaea = len(df_merge_p[df_merge_p["Domain"]=="Archaea"]["Phylum"].unique())
        l_bacteria = len(df_merge_p[df_merge_p["Domain"]=="Bacteria"]["Phylum"].unique())
        len_d_list[0].append(l_archaea)
        len_d_list[1].append(l_bacteria)
        # len_p_list.append(len(dt30[1]))
        name_p_list+=nt30[1]
    # print("len_d_list",len_d_list)
    uni_name_p_list,cnt_uni = np.unique(name_p_list,return_counts=True)
    cnt_uni_std = [ cu/sum(cnt_uni) for cu in cnt_uni ]

    for a_sim in range(cnt_sample):
        # l_ar = random.choice(len_d_list[0])
        # l_bac = random.choice(len_d_list[1])
        
        # len_p = random.choice(len_p_list)
        # len_p = l_ar+l_bac
        
            
        # list_p = np.random.choice(uni_name_p_list,size=len_p,p=cnt_uni_std,replace=False)
        data_sim,dis_sp_sim,dis_shift_sim =get_sp_cnt_data_level_sample(len_d_list,
                            data_tree_30,dis_tree_30,shift_tree_30,position_phy,                                                           uni_name_phy_bac,cnt_name_phy_std_bac,uni_name_phy_ar,cnt_name_phy_std_ar)
        data_sim_root,dis_sim_root,dis_shift_sim_root = get_root_node(data_sim,dis_sp_sim,dis_shift_sim)  
        flag,_ ,_= is_subtree2(data_sim_root,data_tree_root)
        if flag:
            dis_12,name_sp_12 = get_tree_dis_new_3(data_sim_root,dis_sim_root,
                    data_tree_root,dis_tree_root,datad_name_root,0,False,flag_dis)
            name_sp_list.append(list(name_sp_12[-1]))
        
            print(i,"sp",len(name_sp_12[-1]),"dis",dis_12)

else:#一个环境只有一个样本
    
    # for list_sp_genus_ori,data_tree in zip(list_sp_genus_ori_all,data_tree_30):
    for i,data_tree in enumerate(data_tree_30):
        # print("i",i)
        for a_sim in range(cnt_sample):
#             if flag_len_sp:
#                 list_sp_genus_ori = list_sp_genus_ori_all[i]
#                 len_sp_sample = sum(data_tree[-1])
#                 len_sd = 0.1*len_sp_sample
#                 len_sp_new_list = list(range(round(len_sp_sample-len_sd),round(len_sp_sample+len_sd)))
#                 len_sp_new = random.choice(len_sp_new_list)

#                 list_sp_genus = random.sample(list_sp_genus_ori,len_sp_new)
#                 merge_name_sp = taxonomy[taxonomy["Species"].isin(list_sp_genus)].reset_index(drop=True)
#                 data_sim,name_sim = get_num_new_2(merge_name_sp)
#                 dict_sim,data_sim2 = get_dis_list_alltree(merge_name_sp,
#                     path_tax,path_nw,path_sim,path_abs) 
#                 dis_sp_sim,dis_shift_sim = get_dis_tree(dict_sim,name_sim,data_sim,merge_name_sp)
#                 shuffle = False
#             else:
            # print("单环境单样本")
            data_sim = copy.deepcopy(data_tree)#不改物种
            dis_sp_sim = dis_tree_30[i]
            dis_shift_sim = shift_tree_30[i]
            # if cnt_sample==1:
            #     shuffle = False
            # else:
            #     shuffle = True
            index_r0 = 0
            for r0 in data_sim[0]:
                if r0==0:
                    del data_sim[0][index_r0]
            

            data_sim_root,dis_sim_root,dis_shift_sim_root = get_root_node(data_sim,dis_sp_sim,dis_shift_sim) 
            flag,_ ,_= is_subtree2(data_sim_root,data_tree_root)
            
            dis_12,name_sp_12 = get_tree_dis_new_3(data_sim_root,dis_sim_root,
                data_tree_root,dis_tree_root,datad_name_root,0,False,flag_dis)
            name_sp_list.append(name_sp_12[-1])
            print(i,"dis",dis_12)

if path_sim.endswith("/"):
    path_sim_root = path_sim.replace("community/","")
else:
    path_sim_root = path_sim.replace("community","")
path_save = os.path.join(path_sim_root,"species","sp_list.pkl")
print(path_save)
with open(path_save,"wb") as f:
    pickle.dump(name_sp_list,f)
    




