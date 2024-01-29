# -*- coding: utf-8 -*-
# Run Kraken in batch.
# Identify the species in the sample from reads.
import pandas as pd
from pandas import read_csv
import os
import subprocess
import sys



in_file = sys.argv[1]
out_file_root = sys.argv[2]
mode = int(sys.argv[3]) 
dir_kraken = sys.argv[4]
path_db = sys.argv[5]
start = int(sys.argv[6])
end = int(sys.argv[7])
suf = sys.argv[8]
suf_pre = sys.argv[9]


# single-read
def get_report_1(in_file,out_file,path_db,start,end,suffix=".fastq.gz"):
    # single
    
    list_srr = sorted(os.listdir(in_file))[start:end]
    for f in list_srr:
       
        srr = f.replace(suffix,"")
        res=os.path.join(out_file,srr+".out")
        report=os.path.join(out_file,srr+"_rep.txt")
        path_f = os.path.join(in_file,f)
        command = "kraken2 --db "+ path_db+" --memory-mapping "+path_f+" --use-names"+" --output "+res+" --report "+report + " --threads 32 --use-mpa-style"
        
        subprocess.run(command,
                       shell=True,
                       cwd = dir_kraken)


#paired-read
def get_report_2(in_file,out_file,path_db,start,end,suffix1="_1",suffix=".fastq.gz"):
    #paired-end
    list_srr_all = sorted(os.listdir(in_file))
    
    list_srr = [srr.replace(suffix1+suffix,"") for srr in list_srr_all if suffix1 in srr]
    
    for f in list_srr[start:end]:#batch
        
        srr = f.replace(suffix,"")
        
        res=os.path.join(out_file,srr+".out")
        report=os.path.join(out_file,srr+"_rep.txt")
        f1 = f+suffix1+suffix
        f2 = f+suffix1+suffix
        path_f1 = os.path.join(in_file,f1)
        path_f2 = os.path.join(in_file,f2)
        
        
        # To prevent incomplete paired-end data, perform a check.
        a = os.path.exists(path_f1)
        b = os.path.exists(path_f2)

        
        command = "kraken2 --db "+ path_db+" --memory-mapping --paired "+path_f1+" "+path_f2 \
                +" --use-names"+" --output " +res+" --report "+report +" --threads 32 --use-mpa-style"
        
        subprocess.run(command,
                       shell=True,
                       cwd = dir_kraken)
    

out_file = os.path.join(out_file_root,"kraken_res")

    
if mode==1:
    suffix = "."+suf
    get_report_1(in_file,out_file,path_db,start,end,suffix)
else:
    suffix = "."+suf
    get_report_2(in_file,out_file,path_db,start,end,suf_pre,suffix)

if mode == 1:
    srr_batch = sorted(os.listdir(in_file))[start:end]
    
    srr_batch_new = [srr.replace(suffix,"") for srr in srr_batch]
    
else:
    list_srr_all = sorted(os.listdir(in_file))
    
    srr_batch = [srr.replace(suf_pre+suffix,"") for srr in list_srr_all if suf_pre in srr]
    
    srr_batch_new = srr_batch[start:end]

df_batch = pd.DataFrame(srr_batch_new)

path_batch_root = os.path.join(out_file_root,"batch")

path_batch = os.path.join(path_batch_root,str(start)+"_"+str(end)+".txt")
df_batch.to_csv(path_batch,header=None,index = False)

cmd_get_read_cnt = "bash 1_2_read_cnt_abu_newout.sh "+ out_file_root+" "+path_batch
print(cmd_get_read_cnt)

abs_path = os.path.abspath(sys.argv[0])
abs_path_dir = os.path.dirname(abs_path)

subprocess.run(cmd_get_read_cnt,
                shell=True,
                cwd = abs_path_dir)
