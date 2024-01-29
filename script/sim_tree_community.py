# -*- coding: utf-8 -*-
# Generate sequencing data.
import numpy as np
import pandas as pd
import subprocess
import os
import random
import pickle 
import random
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import sys
from scipy.stats import lognorm


path_abs = sys.argv[1]
path_sim = sys.argv[2]
cnt_strain = int(sys.argv[3])
path_abu = sys.argv[4] 
sim_mode = int(sys.argv[5])
genome_strain = int(sys.argv[6]) 

# sim_mode is not distinguished based on env, sample, or community
# 0 indicates that the input abundance is a list, and sampling is done within that list.
#1 indicates that the input is a distribution, used for generation. 
# Sample has two modes: 0 represents multiple samples in one environment, and 1 represents one sample in one environment.
# 2 is a list dictionary. For community, there are two modes: 0 indicates randomly selecting an environment, and 2 indicates that the user has specified the abundance of each species.

if genome_strain == 0:
    path_strain_genome = os.path.join(path_abs,"data/strain_genome/")
else:
    path_strain_genome = os.path.join(path_abs,"data/strain_download/")

def get_strain_abu(depth_sample,list_strain_choice,list_cnt_choice):
    depth_list_strain = [ float("{:.2f}".format(depth_sample[ds_id]/list_cnt_choice[ds_id])) \
     for ds_id in range(len(depth_sample)) for cti in range(list_cnt_choice[ds_id])]
    return depth_list_strain

def sim_strains(genome,depth,sim_out,path_abs,sample_n=1,path_strain_genome = path_strain_genome):
    # Simulate sequencing data using PBSim3.
    # two: method: qshmm and errhmm.
    # path_strain_genome:data/strain_download/
    path_s = os.path.join(sim_out,"strains.csv")
    print(path_s)
    df_g = pd.DataFrame(genome)
    df_a = pd.DataFrame(depth)
    df_s = pd.concat([df_g,df_a],axis = 1)
    df_s = df_s.reset_index(drop=True)
    df_s.columns = ["strain","abu"]
    df_s.to_csv(path_s,sep="\t",index=None)

        
if genome_strain == 0:
    path_sp = os.path.join(path_abs,"data/sp_path.txt")
    df_path_sp = pd.read_csv(path_sp,header=None,sep="\t")
    df_path_sp.columns=["sp","path"]
else:
    path_down = os.path.join(path_abs,"data/strain_file/down_strain.csv")
    df_down = pd.read_csv(path_down,header=None,sep="\t")
    df_down.columns = ["sp","acns","path","genome"]
    df_path_sp = df_down[["sp","genome"]]
    df_path_sp.columns=["sp","path"]


path_sp_sim = os.path.join(path_sim,"species/sp_list.pkl")
path_reads = os.path.join(path_sim)
sim_out = path_reads
print(sim_out)

with open(path_sp_sim,"rb") as f:
    name_sp_list = pickle.load(f)
# List for storing genomes



    


with open(path_abu,"rb") as f:
    abu_real = pickle.load(f)

for n_index,name_sp_sample in enumerate(name_sp_list):# name_sp_sample是一个样本
    sp_a_sample = name_sp_sample[:]
    
    df_sp_sample = pd.DataFrame(sp_a_sample,columns=["sp"])
    merge_sp_path = pd.merge(df_sp_sample,df_path_sp,on=["sp"])
    
    list_cnt_choice = [] 
    # Record the actual number of generated strains for each species, which can be used to control the allocation of abundance (within-species abundance).
    list_strain_choice =[]
    for sp in sp_a_sample:
        
        path_this_sp = merge_sp_path[merge_sp_path["sp"]==sp]["path"].to_list()
        cnt_chioce = min(cnt_strain,len(path_this_sp)) 
        #The number of selected strains must not exceed the maximum number of strains for this species.
        
        list_cnt_choice.append(cnt_chioce)
        strian_choice = random.sample(path_this_sp,cnt_chioce)
        
        strian_choice_new = [os.path.join(path_strain_genome,sc) for sc in strian_choice]
        list_strain_choice += strian_choice_new 
        # This list contains the selected genomes.
    
    len_sp = len(sp_a_sample)

    if sim_mode == 0: # abu_real is a list.
        print("mode 0")
        if len(abu_real)>=1000:
            len_sam = 1000
        else:
            len_sam = len(abu_real)
        abu_sample_sim_1000 = random.sample(abu_real,len_sam)
    elif sim_mode==1: # One simulation corresponds to one abundance.
        print("mode 1")
        abu_log_norm = abu_real[n_index]
        shape_ln = abu_log_norm['shape']
        loc_ln = abu_log_norm['loc']
        scale_ln = abu_log_norm['scale']

    
        abu_sample_sim_1000 = lognorm.rvs(s=shape_ln, loc=loc_ln, scale=scale_ln, size=1000)
    else:
        print("mode 2")
        abu_sample_sim = abu_real[n_index]
    if sim_mode!=2:
        abu_sample_sim = random.sample(abu_sample_sim_1000,len_sp)
        sum_abu_sample =  sum(abu_sample_sim) # Restrict sampling to around 1.
        cnt_while = 0
        flag_re_s = False
        min_sum = 0.98
        max_sum = 1.02
        flag_try_1 = True # Multiple attempts.
        while sum_abu_sample<min_sum or sum_abu_sample>max_sum:
            abu_sample_sim = random.sample(abu_sample_sim_1000,len_sp)
            sum_abu_sample =  sum(abu_sample_sim)
            if cnt_while>100:
                if flag_try_1:
                    flag_try_1 = False # The first attempt failed.
                    min_sum = 0.95
                    max_sum = 1.05
                    cnt_while = 0
                else:

                    flag_re_s = True
                    break
            cnt_while+=1
        if flag_re_s:
            
            abu_sample_sim = random.sample(abu_sample_sim,len_sp)
        
        abu_sample_sim_new = [float("{:.4f}".format(ass)) for ass in abu_sample_sim]
        abu_sample_sim = abu_sample_sim_new[:]
        
        
    # Upon completion of the above execution, the abundance of this sample is obtained.
    print("abu_sample_sim",abu_sample_sim)
    if genome_strain == 1:
        path_tmp_strain = os.path.join(sim_out,"sample_strain_path_"+str(n_index)+".txt")
        
        df_tmp_strain = pd.DataFrame(list_strain_choice)
        df_tmp_strain.columns = ["acns"]
        df_tmp_strain["acns"] = df_tmp_strain["acns"].apply(lambda x:x.split("/")[-1].split("_")[0]+"_"+x.split("/")[-1].split("_")[1])
        df_tmp_strain.to_csv(path_tmp_strain,header=None,sep="\t",index=False)
        path_down = os.path.join(path_abs,"script/down_strains.py")
        cmd_down = "python "+path_down+" "+ path_tmp_strain+" "+path_abs
       
        print("start to download strain genomes")
        try:
            result = subprocess.run(cmd_down,shell=True,cwd="./")
            result.check_returncode()
        except subprocess.CalledProcessError as e:
            sys.exit(1)
        

    min_abu_sample = min(abu_sample_sim)

    # depth_sample = [float("{:.2f}".format(min_depth*(mas/min_abu_sample))) for mas in abu_sample_sim ]\
    depth_sample = abu_sample_sim[:]
    depth_list_strain = get_strain_abu(depth_sample,list_strain_choice,list_cnt_choice)

    print(depth_sample,list_strain_choice,list_cnt_choice)
    print(len(depth_list_strain),depth_list_strain)

    sim_strains(list_strain_choice,depth_list_strain,sim_out,path_abs,n_index+1,path_strain_genome)
        
    
    


