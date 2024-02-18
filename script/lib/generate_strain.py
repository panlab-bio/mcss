# Generate sequence
import pandas as pd
import numpy as np
import os
import random
from Bio import SeqIO
from Bio.Seq import Seq
def get_strain_genome(seq_id_init,record,record_seq_list,ani,name_j,strain_s,path_strian):
    print(ani,name_j,strain_s,strain_s,path_strian)
    len_genome = len(record_seq_list)
    len_base = round(len(record_seq_list)*(1-ani))
    
    
    len_1 = random.randint(round(len_base*0.2),round(len_base*0.8)) 
    len_2 = len_base - len_1 
    
    base_list = ["A","T","C","G"]
    index_variant = random.sample(range(len(record_seq_list)),len_base)
    insert_base=[]
    insert_index=[] 
    for i in range(len(index_variant)):
        if i <len_1:
            
            base_new = random.choice([base for base in base_list if base!=record_seq_list[index_variant[i]]])
            record_seq_list[index_variant[i]] = base_new 
        else:
            flag_indel = random.choice([0,1])
            if flag_indel==0: 
                record_seq_list[index_variant[i]]="N" 
                len_genome-=1
            else: 
                insert_index.append(index_variant[i])
                len_genome+=1
    insert_index.sort(reverse=True)
    
    for i in insert_index:
        base = random.choice(base_list)
        
        record_seq_list.insert(i,base)
    record_seq_list_new = [base for base in record_seq_list if base!="N"]
    separator = ''
    record_seq_list_new_str = separator.join(record_seq_list_new)
    record_seq_list_seq = Seq(record_seq_list_new_str)
    record.seq = record_seq_list_seq
    
    output_file_name = os.path.join(path_strian,name_j)
    
    record.id = seq_id_init+"_"+str(strain_s)
    print("----",strain_s,record.id,output_file_name)
    with open(output_file_name, "w") as output_file:
        SeqIO.write(record, output_file, "fasta")
    
    return record_seq_list_new,len_genome
