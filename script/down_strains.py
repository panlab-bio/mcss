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

path_df = sys.argv[1]
path_abs = sys.argv[2]
df_strain = pd.read_csv(path_df,header=None,sep="\t")
df_strain.columns = ["acns"]


path_down = os.path.join(path_abs,"data/strain_file/down_strain.csv")
# print("path_down\n",path_down)
df_down = pd.read_csv(path_down,header=None,sep="\t")
df_down.columns = ["sp","acns","path","genome"]
df_merge = pd.merge(df_strain,df_down,on=["acns"])
# print("df_merge\n",df_merge)

path_strain_genome = os.path.join(path_abs,"data/strain_download")
list_saved = os.listdir(path_strain_genome)
# print("list_saved\n",list_saved)
mask = ~df_merge["genome"].isin(list_saved)
df_mask = df_merge[mask]
list_path = df_mask["path"].dropna().tolist()
# print(list_path)


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
        
# path_cur = os.path.dirname(path_abs)
# path_cur = "../data/"
path_strain_genome_dl = os.path.join(path_abs,"data/strain_genome_dl")
os.makedirs(path_strain_genome_dl,exist_ok=True)


    
for lp in list_path:
    # print(lp)
    print("download ",lp)
    response = requests.get(lp)
    if response.status_code == 200:
        file_name = os.path.join(path_strain_genome_dl, os.path.basename(lp))
        with open(file_name, 'wb') as file:
            file.write(response.content)
            
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
