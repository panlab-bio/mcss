import requests
import os
import subprocess
import pandas as pd
import numpy as np
import os
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import random
import sys
path_abs = os.path.abspath(sys.argv[0])

path_cur = os.path.dirname(path_abs)
path_mcss = os.path.dirname(path_cur)
print(path_cur)
print(path_mcss)
def replace_n_with_random(sequence):
    replaced_sequence = ""
    nucleotides = ['A', 'T', 'C', 'G']
    for base in sequence:
        if base not in nucleotides:
            replaced_sequence += random.choice(nucleotides)
        else:
            replaced_sequence += base
    return replaced_sequence

def process_fasta_file(input_file, output_file_name):
    with gzip.open(input_file, "rt") as fasta_gz_file:
        # has_n = any("N" in record.seq for record in SeqIO.parse(fasta_gz_file, "fasta"))
        # print("has_n",has_n)
        records = []
        # if has_n:
        cnt=0
        for record in SeqIO.parse(fasta_gz_file, "fasta"):
            if "N" in record.seq:
                replaced_sequence = replace_n_with_random(str(record.seq))
                new_record = record
                new_record.seq = Seq(replaced_sequence)
                records.append(new_record)
                # print("N",cnt)
            else:
                records.append(record)
                # print("no N",cnt)
            cnt+=1

        # print("records",len(records))

        with open(output_file_name, "w") as output_file:
            SeqIO.write(records, output_file, "fasta")
    
#把结果放到一起
def concatenate_sequences(input_file, output_file_name):
    sequences = []
    flag = True
    
    for record in SeqIO.parse(input_file, "fasta"):
        if flag:
            record_1 = record
            flag = False
        sequences.append(str(record.seq))
    concatenated_sequence = ''.join(sequences)
    record_1.seq=Seq(concatenated_sequence)
    with open(output_file_name,"w") as output_file:
        SeqIO.write(record_1,output_file,"fasta")

def acns_path_fun_new_2(acn_a,path_genome):
# acn_a = "RS_GCF_004114995.1"|
    name_acn = acn_a.split("_")
    # print(name_acn)
    path_a = os.path.join(path_genome,name_acn[0],name_acn[1][:3],name_acn[1][3:6],name_acn[1][6:9])
    genome = os.path.join(path_a,os.listdir(path_a)[0])
    return genome



path_gt = "https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_reps/gtdb_genomes_reps_r207.tar.gz"
cmd_gt = "wget -c "+path_gt
subprocess.run(cmd_gt,shell=True,cwd="./")
cmt_tar = "tar -xzf ./gtdb_genomes_reps_r207.tar.gz"
subprocess.run(cmt_tar,shell=True,cwd="./")

path_genome = "./gtdb_genomes_reps_r207"


path_strain_genome = os.path.join(path_mcss,"data/strain_genome")
os.makedirs(path_strain_genome,exist_ok=True)
path_genome = "/data/huixingqi/sim_meta/data/ncbi_genome/gtdb_genomes_reps_r207/"
path_acns_gt = os.path.join(path_mcss,"data/gtdb_taxonomy_acns.tsv")
df_acns_gt = pd.read_csv(path_acns_gt,header=None,sep="\t")
df_acns_gt.columns=["acns","tax"]
df_acns_gt["acns"] = df_acns_gt["acns"].apply(lambda x:x.replace("GB_",""))
df_acns_gt["acns"] = df_acns_gt["acns"].apply(lambda x:x.replace("RS_",""))
list_gt_acns = df_acns_gt["acns"].to_list()
print(len(list_gt_acns))
# for lga in list_gt_acns[:2]: 
for lga in list_gt_acns: 
    acns_path = acns_path_fun_new_2(lga,path_genome)
    
    print(acns_path)
    strain_name = acns_path.split("/")[-1].replace(".gz","")
    path_strain_name = os.path.join(path_strain_genome,strain_name)
    if not os.path.exists(path_strain_name):
        process_fasta_file(acns_path,path_strain_name)
        concatenate_sequences(path_strain_name,path_strain_name)




path_strain_path = os.path.join(path_mcss,"data/strain_file/strain_genome_path.txt")

df_strain_path = pd.read_csv(path_strain_path,header=None,sep="\t")
df_strain_path.columns= ["path"]

path_strain_genome_dl = os.path.join(path_cur,"strain_genome_dl")
os.makedirs(path_strain_genome_dl,exist_ok=True)

list_path = df_strain_path["path"].to_list()




# for lp in list_path[:2]:
for lp in list_path:
    print(lp)
    response = requests.get(lp)
    if response.status_code == 200:
        file_name = os.path.join(path_strain_genome_dl, os.path.basename(lp))
        with open(file_name, 'wb') as file:
            file.write(response.content)
            
#这个基因组用来保存解压后的fasta文件

list_strain_choice_all = os.listdir(path_strain_genome_dl)
                
for strain in list_strain_choice_all:
    # print(strain)
    
    strain_name = strain.split("/")[-1].replace(".gz","")
    # print(strain_name)
    strain_new = os.path.join(path_strain_genome_dl,strain)
    
    path_strain_name = os.path.join(path_strain_genome,strain_name)
    # print(strain_new,path_strain_name)
    if not os.path.exists(path_strain_name):
    #     # print(path_strain_name)
        # cmd_unzip = "gunzip -c "+ strain +" > "+ path_strain_name
        # subprocess.run(cmd_unzip,shell=True)
        #根本不需要解压，因为拼接和去N操作输入就是gz,输出就是非压缩文件
        
        process_fasta_file(strain_new,path_strain_name)
        concatenate_sequences(path_strain_name,path_strain_name)
        
            
path_sp_strain_path = os.path.join(path_mcss,"data/strain_file/sp_strain.txt")
path_sp_path =  os.path.join(path_mcss,"data/sp_path.txt")

list_all = os.listdir(path_strain_genome)
#all 
df_list_all = pd.DataFrame(list_all,columns=["path"])
df_list_all["path"] = df_list_all["path"].apply(lambda x:x.replace("_genomic.fna",""))

#download
df_list_dl = pd.read_csv(path_sp_strain_path,header=None,sep="\t")
df_list_dl.columns=["sp","path"]
df_list_dl["path"] = df_list_dl["path"].apply(lambda x:x.split("/")[-1])

df_list_gt = pd.read_csv(path_acns_gt,header=None,sep="\t")
df_list_gt.columns=["path","tax"]
df_list_gt["path"] = df_list_gt["path"].apply(lambda x:x.replace("GB_",""))
df_list_gt["path"] = df_list_gt["path"].apply(lambda x:x.replace("RS_",""))
df_list_gt["sp"] = df_list_gt["tax"].apply(lambda x:x.split(";")[-1].replace(" ","_"))

df_list_gt_sp = df_list_gt[["sp","path"]]
df_list_dl_sp = df_list_dl[["sp","path"]]
df_list_gt_dl = pd.concat([df_list_gt_sp,df_list_dl_sp],axis=0,ignore_index=True)

df_sp_all_merge = pd.merge(df_list_gt_dl,df_list_all,on=["path"])
df_sp_all_merge["path"] = df_sp_all_merge["path"].apply(lambda x:x+"_genomic.fna")
print(len(df_sp_all_merge))
df_sp_all_merge.to_csv(path_sp_path,header=None,sep="\t",index=None)



