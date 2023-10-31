# -*- coding: utf-8 -*-
#批量跑kraken
# nohup python 1_kraken.py marine 0 10 2&>log/run_1_kraken.log &

import pandas as pd
from pandas import read_csv
import os
import subprocess
import sys

#首先要判断是什么模式，有的是双端，有的是单端
#需要输入文件夹和输出文件夹,输出文件根据输入生成，不用指定了
# env = sys.argv[1] #
#这个版本就不同环境目录不一样了
in_file = sys.argv[1]
out_file_root = sys.argv[2]
# mode = 2 #默认是双端哈，不默认了，反正也不能省略
mode = int(sys.argv[3]) #是什么模式，mode==1是单端，2是双端
dir_kraken = sys.argv[4]
path_db = sys.argv[5]
start = int(sys.argv[6])
end = int(sys.argv[7])
# path_db = "/dev/shm/db/"
# in_file = os.path.join("/data/huixingqi/data/env",env)


def get_report_1(in_file,out_file,path_db=path_db):
    #单端
    # path_db = "/dev/shm/db/"
    list_srr = sorted(os.listdir(in_file))
    for f in list_srr:
        srr = f.split(".")[0]
        res=os.path.join(out_file,srr+".out")
        report=os.path.join(out_file,srr+"_rep.txt")
        path_f = os.path.join(in_file,f)
        command = "kraken2 --db "+ path_db+" --memory-mapping "+path_f+" --use-names"+" --output "+res+" --report "+report + " --threads 32 --use-mpa-style"
        # print(command)
        subprocess.run(command,
                       shell=True,
                       cwd = dir_kraken)



def get_report_2(in_file,out_file,path_db=path_db):

    list_srr_all = sorted(os.listdir(in_file))
    #双端的问题是有两个srr，_1和_2，所以要获得根name
    list_srr = [srr.replace("_1.fastq.gz","") for srr in list_srr_all if "_1" in srr]
    suffix = list_srr_all[0].split("_")[-1][1:]
    # print(suffix)
    # print(list_srr)
    for f in list_srr:#batch
        # srr = f.split("_")[0]
        srr = f.split(".")[0]
        # print(srr)
        res=os.path.join(out_file,srr+".out")
        report=os.path.join(out_file,srr+"_rep.txt")
        f1 = f+"_1"+suffix
        f2 = f+"_2"+suffix
        path_f1 = os.path.join(in_file,f1)
        path_f2 = os.path.join(in_file,f2)
        # print(path_f1,path_f2)
        
        #为了防止有的双端数据不全，判断一下
        a = os.path.exists(path_f1)
        b = os.path.exists(path_f2)

        # print(res)
        command = "kraken2 --db "+ path_db+" --memory-mapping --paired "+path_f1+" "+path_f2 \
                +" --use-names"+" --output " +res+" --report "+report +" --threads 32 --use-mpa-style"
        # print(command)
        subprocess.run(command,
                       shell=True,
                       cwd = dir_kraken)
    

out_file = os.path.join(out_file_root,"kraken_res")
#放在py中生成批量时会出错，所以放在外面跑
# if os.path.exists(out_file)==False:
#     os.makedirs(out_file)
# print(in_file)
# print(out_file)
    
if mode==1:
     get_report_1(in_file,out_file)
else:
     get_report_2(in_file,out_file)

if mode == 1:
    srr_batch = sorted(os.listdir(in_file))[start:end]
    srr_batch_new = [srr.split(".")[0] for srr in srr_batch]
else:
    list_srr_all = sorted(os.listdir(in_file))
    srr_batch = [srr.split(".")[0].replace("_1","") for srr in list_srr_all if "_1" in srr]
    srr_batch_new = srr_batch[start:end]

df_batch = pd.DataFrame(srr_batch_new)
# path_batch_root = os.path.join("/data/huixingqi/data/env/res/batch",env)
path_batch_root = os.path.join(out_file_root,"batch")
# if not os.path.exists(path_batch_root):
#     os.makedirs(path_batch_root)
path_batch = os.path.join(path_batch_root,str(start)+"_"+str(end)+".txt")
df_batch.to_csv(path_batch,header=None,index = False)

cmd_get_read_cnt = "bash 1_2_read_cnt_abu_newout.sh "+ out_file_root+" "+path_batch
print(cmd_get_read_cnt)

abs_path = os.path.abspath(sys.argv[0])
abs_path_dir = os.path.dirname(abs_path)
# print(abs_path_dir)
subprocess.run(cmd_get_read_cnt,
                shell=True,
                cwd = abs_path_dir)