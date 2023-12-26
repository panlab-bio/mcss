# -*- coding: utf-8 -*-
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
#这个专门为sample模式下生成测序数据

path_abs = sys.argv[1]
path_sim = sys.argv[2]#这个的输入是真的sim，不是community，之前的因为牵连太多，就用community了
cnt_strain = int(sys.argv[3])
depth = int(sys.argv[4])
type_depth = int(sys.argv[5])
path_abu = sys.argv[6] #丰度分布文件
sim_mode = int(sys.argv[7])
# path_pbsim = sys.argv[8] 
# method = "qshmm"
method = sys.argv[8] 
# model = "/data/huixingqi/software/1_tools/pbsim3/data/QSHMM-RSII.model"
model =  sys.argv[9]
pass_sum = sys.argv[10]
genome_strain = int(sys.argv[11]) # 0 or 1,0是早就存好了，1是现下

#sim_mode不是按照env,sample,community来区分
#0代表输入的abu是一个列表，在里面采样
#1代表输入的是分布，用于generate 。sample有两种模式 0 代表一个环境多个样本，1代表一个环境一个样本
#2是一个列表字典，community有两种模式 0 是值得随机找个环境的用 2表示用户指定了每个物种的丰度

# print("pbsim-----------",path_pbsim)
#type_depth为0的时候，认为输入的是最小depth
#否则认为输入的是平均depth
if type_depth==0:
    print("min")
    flag_min_depth = True
    min_depth = depth
else:
    print("mean")
    flag_min_depth = False
    mean_depth = depth



def get_strain_abu(depth_sample,list_strain_choice,list_cnt_choice):
    depth_list_strain = [ float("{:.2f}".format(depth_sample[ds_id]/list_cnt_choice[ds_id])) \
     for ds_id in range(len(depth_sample)) for cti in range(list_cnt_choice[ds_id])]
    return depth_list_strain

def sim_read_pbsim3(genome,depth,method,model,sim_out,path_abs,sample_n=1,l_min=100,l_max=1000000):
    #用pbsim3模拟测序数据
    # method有两种，qshmm 和 errhmm
    path_concat = os.path.join(sim_out,"sim_concat/sim.fastq")
    
    if os.path.exists(path_concat):
        print(path_concat)
        subprocess.run("rm "+path_concat,shell=True)
    # cmd_pbsim_head = "bash run_pbsim3.sh "+sim_out+" "+method+" "+model+" "+path_pbsim
    cmd_pbsim_head = "bash run_pbsim3.sh "+sim_out+" "+method+" "+model #因为pbsim内置到环境中了，所以不需要指定pbsim路径了
    # print(cmd_pbsim_head)
    my_pre_head="sim_"
    cnt_sim=1
    for my_g,my_d in zip(genome,depth):
        my_pre=my_pre_head+str(cnt_sim)
        cnt_sim+=1
        cmd_pbsim = cmd_pbsim_head+ " "+str(my_d)+" "+my_g+" "+my_pre+" "+pass_sum+" "+str(sample_n)
        path_shell=os.path.join(path_abs,"script")
        # print(path_shell)
        # print(cmd_pbsim)
        subprocess.run(cmd_pbsim,shell=True,cwd=path_shell)
        # try:
        #     result = subprocess.run(cmd_pbsim, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,cwd="./")
        #     print("Output:", result.stdout)
        # except subprocess.CalledProcessError as e:
        #     print("Error:", e)
        #     print("Output:", e.stdout)
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
# print(df_path_sp["sp"].nunique())

path_sp_sim = os.path.join(path_sim,"species/sp_list.pkl")
path_reads = os.path.join(path_sim,"pbsim")
sim_out = path_reads
print(sim_out)

with open(path_sp_sim,"rb") as f:
    name_sp_list = pickle.load(f)
# print(name_sp_list)
#存放基因组的列表


if genome_strain == 0:
    path_strain_genome = os.path.join(path_abs,"data/strain_genome/")
else:
    path_strain_genome = os.path.join(path_abs,"data/strain_download/")
    

    
    
    
    
    
# if sim_mode!=2:
with open(path_abu,"rb") as f:
    abu_real = pickle.load(f)
# else:
#     abu_real = pd.read_csv(path_abu,header=None,sep="\t")
#     abu_real.columns=["abu"]
#     abu_sample_sim = abu_real["abu"].to_list()
for n_index,name_sp_sample in enumerate(name_sp_list):# name_sp_sample是一个样本
    sp_a_sample = name_sp_sample[:]
    # print(sp_a_sample)
    df_sp_sample = pd.DataFrame(sp_a_sample,columns=["sp"])
    merge_sp_path = pd.merge(df_sp_sample,df_path_sp,on=["sp"])
    
    list_cnt_choice = [] #记录每个物种的实际生成的strain个数，这样可以用于控制丰度的分配（物种内丰度）
    list_strain_choice =[]
    for sp in sp_a_sample:
        # print(sp)
        path_this_sp = merge_sp_path[merge_sp_path["sp"]==sp]["path"].to_list()
        cnt_chioce = min(cnt_strain,len(path_this_sp)) #选的strain个数不能超过这个物种的strain的个数的最大值
        # print("cnt_chioce",cnt_chioce)
        list_cnt_choice.append(cnt_chioce)
        strian_choice = random.sample(path_this_sp,cnt_chioce)
        # down_strains(strian_choice)
        strian_choice_new = [os.path.join(path_strain_genome,sc) for sc in strian_choice]
        list_strain_choice += strian_choice_new #这个列表里面存的是选好的基因组
    # list_strain_choice是选好的基因组
    len_sp = len(sp_a_sample)

    if sim_mode == 0:#这时候abu_real是一个大列表
        print("mode 0")
        if len(abu_real)>=1000:
            len_sam = 1000
        else:
            len_sam = len(abu_real)
        abu_sample_sim_1000 = random.sample(abu_real,len_sam)
    elif sim_mode==1:#这时候一次模拟，对应一个丰度
        print("mode 1")
        abu_log_norm = abu_real[n_index]
        shape_ln = abu_log_norm['shape']
        loc_ln = abu_log_norm['loc']
        scale_ln = abu_log_norm['scale']

        # dist = abu_real[n_index]
        # abu_sample_sim_1000 = list(dist.generate(1000))
        abu_sample_sim_1000 = lognorm.rvs(s=shape_ln, loc=loc_ln, scale=scale_ln, size=1000)
    else:
        print("mode 2")
        abu_sample_sim = abu_real[n_index]
    if sim_mode!=2:
        abu_sample_sim = random.sample(abu_sample_sim_1000,len_sp)
        sum_abu_sample =  sum(abu_sample_sim)#限制采样的和在1附近
        cnt_while = 0
        flag_re_s = False
        min_sum = 0.98
        max_sum = 1.02
        flag_try_1 = True #表示第一次尝试，第一次尝试不行，第二次尝试，这时候把阈值改一下
        while sum_abu_sample<min_sum or sum_abu_sample>max_sum:
            abu_sample_sim = random.sample(abu_sample_sim_1000,len_sp)
            sum_abu_sample =  sum(abu_sample_sim)
            if cnt_while>100:
                if flag_try_1:
                    flag_try_1 = False #第一次尝试失败
                    min_sum = 0.95
                    max_sum = 1.05
                    cnt_while = 0
                else:

                    flag_re_s = True
                    break
            cnt_while+=1
        if flag_re_s:
            # print("没找到合适的，随便找吧")
            abu_sample_sim = random.sample(abu_sample_sim,len_sp)
        # elif flag_try_1:
        #     # print("第一次找到合适的了",sum_abu_sample)
        # else:
            # print("第二次找到合适的了",sum_abu_sample)
        abu_sample_sim_new = [float("{:.4f}".format(ass)) for ass in abu_sample_sim]
        abu_sample_sim = abu_sample_sim_new[:]
        # sum_abu = sum(abu_sample_sim_new)
        
    #上面执行完，就得到了这个样本的丰度
    print("abu_sample_sim",abu_sample_sim)
    if genome_strain == 1:
        path_tmp_strain = os.path.join(sim_out,"sample_strain_path_"+str(n_index)+".txt")
        # print(path_tmp_strain)
        df_tmp_strain = pd.DataFrame(list_strain_choice)
        df_tmp_strain.columns = ["acns"]
        df_tmp_strain["acns"] = df_tmp_strain["acns"].apply(lambda x:x.split("/")[-1].split("_")[0]+"_"+x.split("/")[-1].split("_")[1])
        df_tmp_strain.to_csv(path_tmp_strain,header=None,sep="\t",index=False)
        path_down = os.path.join(path_abs,"script/down_strains.py")
        cmd_down = "python "+path_down+" "+ path_tmp_strain+" "+path_abs
        # print(cmd_down)
        print("start to download strain genomes")
        subprocess.run(cmd_down,shell=True,cwd="./")
        
    if flag_min_depth:#最小depth
        print("min depth")
        min_abu_sample = min(abu_sample_sim)
        
        depth_sample = [float("{:.2f}".format(min_depth*(mas/min_abu_sample))) for mas in abu_sample_sim ]
        depth_list_strain = get_strain_abu(depth_sample,list_strain_choice,list_cnt_choice)
        # print(depth_list_strain)
        print(depth_sample,list_strain_choice,list_cnt_choice)
        print(len(depth_list_strain),depth_list_strain)
        # print(len(list_strain_choice),list_strain_choice)
        sim_read_pbsim3(list_strain_choice,depth_list_strain,method,model,sim_out,path_abs,n_index+1)
    else:
        print("mean depth")
        sum_depth = len_sp * mean_depth
        depth_sample = [float("{:.2f}".format(sum_depth*mas)) for mas in abu_sample_sim ]
        depth_list_strain = get_strain_abu(depth_sample,list_strain_choice,list_cnt_choice)
        sim_read_pbsim3(list_strain_choice,depth_list_strain,method,model,sim_out,path_abs,n_index+1)
    
    


