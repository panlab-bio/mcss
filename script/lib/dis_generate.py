import numpy as np
import pandas as pd
import itertools
import subprocess
import os
import random
from scipy import stats
# import pickle 
import copy
import json
# get dis list of microbial community
def get_distance_dict_node_alltree_2(n,merge_name,my_columns,dis_n):
    # Return the distances between layers; the input is a specific layer.
    clm=my_columns[n]
    len_my_columns = len(my_columns)
    # Previously, dis_l stored the distance between each lineage and layer. Now, it has been transformed to store distances between nodes, with the key being the lineage and the value being a list.
    dis_l = dict()
    dis_list=[]
    # As soon as a new entity point appears, add it to this list. When adding a new point, check if the point above it is in this list.
    list_above = [] 
    cnt_n=1
    while (n+cnt_n)<len_my_columns:
        dis_list.append([])
        cnt_n+=1
    
    for i in dis_n[n].keys():#i is a specific lineage, such as the bacterial domain

        list_dis_i = []# Record the node distance list for this lineage.
        list_sp_i_node = [] # Record species with entity points above them
        cmd_lineage="merge_name[merge_name."+my_columns[n]+"=="+"i]"
        # Obtain the dataframe for this lineage.
        lineage = eval(cmd_lineage)
        
        cnt_n = 1 # 1 next layer，2 the layer after the next
        
        # dis_n 0 domain 1 phylum 2 class 3order 4 family 5 genus
        if n==0 and len(lineage)==1:
            list_dis_i = [{},{},{},{},{}]
        else:
            while (n+cnt_n)<len_my_columns-1:
                list_dis_next = [] # Distance between lineage i and a node at a specific layer.
                list_dis_next_dict = dict()
                cmd_genus="lineage."+my_columns[n+cnt_n]+".unique()"

                genus = eval(cmd_genus) # Unique list
                
                if len(genus)>1 or n==0:
                    # entity node
                    
                    for j in genus:# Specific taxonomic units.

                        if j in dis_n[n+cnt_n].keys() :# This is also an entity point.
                            cmd_lineage_j="merge_name[merge_name."+my_columns[n+cnt_n]+"=="+"j]"
                            # The dataframe corresponding to the specific phylum j.
                            lineage_j = eval(cmd_lineage_j)

                            cnt_below = cnt_n+1 # The layer below the next layer (down two layers).
                            while (n+cnt_below)<=len_my_columns-1:

                                cmd_genus_j="lineage_j."+my_columns[n+cnt_below]+".unique()"
                                genus_j = eval(cmd_genus_j) # Unique list
                                cnt_below+=1

                            
                            s1 = list(dis_n[n+cnt_n][j].keys())[0]
                            
                            a_dis = np.float64('%.6f'%(dis_n[n][i][s1]-dis_n[n+cnt_n][j][s1]))

                            list_dis_next.append(a_dis)
                            list_dis_next_dict[j] = a_dis
                            dis_list[cnt_n-1].append(a_dis)

                list_dis_i.append(list_dis_next_dict)
                cnt_n+=1
        dis_i_sp = [ key for key,val in dis_n[n][i].items()]
        list_dis_next_sp = []
        list_dis_next_dict_sp = dict()
        for keysp in dis_i_sp:
            
            a_dis_sp = dis_n[n][i][keysp]
            list_dis_next_sp.append(a_dis_sp)
            list_dis_next_dict_sp[keysp] = a_dis_sp
            dis_list[-1].append(a_dis_sp)
        list_dis_i.append(list_dis_next_dict_sp)
        dis_l[i] = list_dis_i
    return dis_l,dis_list
# dis list of species to the lineage eg: The distance of each species within a phylum to this phylum.
def get_distance2(i,merge_name,my_columns,taxonomy,path_nw,path_sim,abs_path_dir):
    dis = []
    # i refers to the column with index i, where indexing starts from zero.
    cmd_unique = "merge_name."+my_columns[i]+".unique()"
    l_unique = eval(cmd_unique) 
    
    cnt=0
    ss = 0
    for j in l_unique:
        domain = merge_name[merge_name[my_columns[i]]==j][my_columns[0]].unique()[0]
        
        cmd_select = "merge_name[merge_name."+my_columns[i] + "==j]"
        aa = eval(cmd_select) # `exec` does not return values, use `eval` instead.       
        len_s = len(aa) # species number
        cmd_len = "len(aa."+my_columns[i+1]+".unique())"
        len_aa = eval(cmd_len)
        ss+=len_s 
        if len_aa>1 : 
            aa_s = aa.Species
            cnt+=len_s
            cmd_select_tax = "taxonomy[taxonomy."+my_columns[i] + "==j]"
            aa_tax = eval(cmd_select_tax)
            aa_s_tax = aa_tax.Species
            dis_g_tax = get_subtree_new(aa_s_tax,domain,path_nw,path_sim,abs_path_dir)
            dis_g = []
            for aa_s_i in aa_s:
                dis_g.append(dis_g_tax[aa_s_i])
            dis+=dis_g
    
    
    return dis
    
# dis dict of species to the lineage, key: sp ,val : dis
def get_subtree_new(a,domain,path_nw,path_sim,abs_path_dir):
    #return dict, key sp, val dis 
    path_data = os.path.join(abs_path_dir,"data")
    path_tg = os.path.join(path_sim,"tmp_data","tmp_g.txt")
    path_tg_root = os.path.join(path_sim,"tmp_data")
    a.to_csv(path_tg,index=False,header=None) 
    
    if domain=="Bacteria":
        tree_name = os.path.join(path_data,"gtdb_tree/bac_sp_new.tree")
    else:
        tree_name = os.path.join(path_data,"gtdb_tree/ar_sp_new.tree")
        
    
    path_prune = os.path.join(path_nw,"nw_prune")
    cmd_sub = path_prune+" -v "+ tree_name+" -f tmp_g.txt >sub_g.txt"
    
    
    subprocess.run(cmd_sub,shell=True,cwd=path_tg_root)
    
    path_distance = os.path.join(path_nw,"nw_distance")
    cmd_dis = path_distance+" -n sub_g.txt"
    sp = subprocess.run(cmd_dis,shell=True,cwd=path_tg_root,
                        capture_output=True,encoding='utf-8')
    
    dis_g = sp.stdout.split()

    idx=0
    dis_dict=dict()
    for i in range(int(len(dis_g)/2)):
        
        dis_dict[dis_g[idx]]=np.float64(dis_g[idx+1])
        idx+=2             
    return dis_dict
    
# dis dictionary, key: lineage(eg phylum) ,val: dict(key:sp, val:dis)
def get_distance_new(i,merge_name,my_columns,taxonomy,path_nw,path_sim,abs_path_dir):
    # return dict, key taxon (eg phylum) , val dis of sp to taxon 
    dis = dict()
    cmd_unique = "merge_name."+my_columns[i]+".unique()"
    l_unique = eval(cmd_unique)
    
    cnt=0
    ss = 0
    for j in l_unique:
        domain = merge_name[merge_name[my_columns[i]]==j][my_columns[0]].unique()[0]
        
        cmd_select = "merge_name[merge_name."+my_columns[i] + "==j]"
        aa = eval(cmd_select) 
        len_s = len(aa)
        cmd_len = "len(aa."+my_columns[i+1]+".unique())"
        len_aa = eval(cmd_len)
        ss+=len_s
        
        if len_aa>1 :
            aa_s = aa.Species
            cnt+=len_s
            cmd_select_tax = "taxonomy[taxonomy."+my_columns[i] + "==j]"
            aa_tax = eval(cmd_select_tax)
            aa_s_tax = aa_tax.Species
            dis_g_tax = get_subtree_new(aa_s_tax,domain,path_nw,path_sim,abs_path_dir)
            dis_g = dict()
            for aa_s_i in aa_s:
                dis_g[aa_s_i] = dis_g_tax[aa_s_i]
            dis[str(j)] = dis_g
        else:
            if i==0 and len_aa==1:  # domain
                
                dis_g = dict()
                aa_s = aa.Species
                
                cmd_select_tax = "taxonomy[taxonomy."+my_columns[i] + "==j]"
                aa_tax = eval(cmd_select_tax)
                aa_s_tax = aa_tax.Species
                dis_g_tax = get_subtree_new(aa_s_tax,domain,path_nw,path_sim,abs_path_dir)
                for aa_s_i in aa_s:
                    dis_g[aa_s_i] = dis_g_tax[aa_s_i]
                dis[str(j)] = dis_g
            
    
    return dis

# Obtain a list of the number of microbial community 
def get_num_new_2(merge_name):
    
    data = merge_name
    
    clu = data.columns
    datad =[y for y in ([data[data[clu[i]]==x][clu[i+1]].nunique() \
                         for x in data[clu[i]].unique()] for i in range(len(clu)-1))]
    datad_name = [y for y in ([x \
                         for x in data[clu[i]].unique()] for i in range(len(clu)))]
    return datad,datad_name

# gtdb ref dis 
def get_dis_list_alltree(df_merge,path_tax,path_nw,path_sim,abs_path_dir):
    #dis of node to node 

    
    taxonomy =pd.read_csv(path_tax,sep="\t")
    
    my_columns = df_merge.columns
    
    data_sp = []
    dis_n=[]
    for i in range(6):
        
        dis_i = get_distance_new(i,df_merge,my_columns,taxonomy,path_nw,path_sim,abs_path_dir)
        dis_n.append(dis_i) # dis of sp to lineage
    
    dis_5 = get_distance2(5,df_merge,my_columns,taxonomy,path_nw,path_sim,abs_path_dir)
    
    dis_g_new = dis_n[5] 
    dis_g_list =dis_5

    dis_f_new,dis_f_list = get_distance_dict_node_alltree_2(4,df_merge,my_columns,dis_n)
    dis_o_new,dis_o_list = get_distance_dict_node_alltree_2(3,df_merge,my_columns,dis_n)
    dis_c_new,dis_c_list = get_distance_dict_node_alltree_2(2,df_merge,my_columns,dis_n)
    dis_p_new,dis_p_list = get_distance_dict_node_alltree_2(1,df_merge,my_columns,dis_n)
    dis_d_new,dis_d_list = get_distance_dict_node_alltree_2(0,df_merge,my_columns,dis_n)
    dis_all_dict = dict()
    dis_all_dict["d"]=[dis_d_new,dis_d_list]
    dis_all_dict["p"]=[dis_p_new,dis_p_list]
    dis_all_dict["c"]=[dis_c_new,dis_c_list]
    dis_all_dict["o"]=[dis_o_new,dis_o_list]
    dis_all_dict["f"]=[dis_f_new,dis_f_list]
    dis_all_dict["g"]=[dis_g_new,dis_g_list]
    return dis_all_dict,data_sp

# Determine the upper-level entity node of the given node.
def get_up_node_cnt(l1,n,taxonomy):
    clu = taxonomy.columns
    cnt_up = 1
    cnt_level = 0
    name_up = l1
    while n-1>=0 and cnt_up<=1:
        name_up = taxonomy[taxonomy[clu[n]] == name_up][clu[n-1]].unique()[0] # upon lineage
        
        cnt_up = taxonomy[taxonomy[clu[n-1]] == name_up][clu[n]].nunique()
        n-=1
        cnt_level+=1
    
    return cnt_level,name_up

# Retrieve the evolutionary distance list of the multiway tree
def get_dis_tree(dis_all_dict_alltree,datad_name,data_sp_alltree,taxonomy):
    cnt_rank = 0
    dis_tree = []
    shift_tree = []
    
    clu = list(dis_all_dict_alltree.keys())
    for rank in datad_name:
        dis_tree_rank = []
        shift_tree_rank = []
        if cnt_rank==0:

            dis_tree_rank = [0 for r0 in rank]
            shift_tree_rank = [0 for r0 in rank]

        else:
            
            cnt_l = 0
            for l in rank:
                if cnt_rank<len(data_sp_alltree):
                    if data_sp_alltree[cnt_rank][cnt_l]>1:
                        
                        cnt_level,name_up = get_up_node_cnt(l,cnt_rank,taxonomy)
                        
                        dis_l = dis_all_dict_alltree[clu[cnt_rank-cnt_level]][0][name_up][cnt_level-1][l]
                        
                        dis_tree_rank.append(dis_l)
                        shift_tree_rank.append(cnt_level-1) 

                    else:
                        cnt_level,name_up = get_up_node_cnt(l,cnt_rank,taxonomy)
                        dis_tree_rank.append(-1)
                        
                        shift_tree_rank.append(cnt_level-1)
                else:
                    cnt_level,name_up = get_up_node_cnt(l,cnt_rank,taxonomy)
                    
                    if clu[cnt_rank-cnt_level]=="g":#
                        dis_l = dis_all_dict_alltree[clu[cnt_rank-cnt_level]][0][name_up][l]
                    else:
                        dis_l = dis_all_dict_alltree[clu[cnt_rank-cnt_level]][0][name_up][cnt_level-1][l]
                    dis_tree_rank.append(dis_l)
                    shift_tree_rank.append(cnt_level-1)
                cnt_l+=1
        cnt_rank+=1
        dis_tree.append(list(dis_tree_rank))   
        shift_tree.append(list(shift_tree_rank))
    return dis_tree,shift_tree

# add root node
def get_root_node(data_sim,dis_sp_sim,dis_shift_sim):
    
    if len(data_sim[0])==1:
        data_sim_root = [[1]]+data_sim
        dis_sp_sim[0][0] = 0  
    else:
        data_sim_root = [[2]]+data_sim
        dis_sp_sim[0][0] = 0
        dis_sp_sim[0][1] = 0
    
    dis_sp_sim_root = [[0]]+dis_sp_sim
    dis_shift_sim_root = [[0]]+dis_shift_sim
    return data_sim_root,dis_sp_sim_root,dis_shift_sim_root

def get_dis_bottom2root_sim(data_list,dis_list,dis_shift_sim):

    dis_bottom2root = copy.deepcopy(dis_list)
    cnt_rank = 0 
    
    for rank in dis_bottom2root[:-1]:
        cnt_l = 0 
        if cnt_rank>=1: 
            
            for rank_l in rank:
                
                if rank_l>0: 

                    start = sum(data_list[cnt_rank][:cnt_l]) 
                    
                    end = start+data_list[cnt_rank][cnt_l]
                    
                    sum_end = 0
                    for se in range(start,end): 
            
                        if dis_bottom2root[cnt_rank+1][se]>0:
           
                            dis_bottom2root[cnt_rank+1][se]+=rank_l 
                        if cnt_rank+1<len(data_list):
                            sum_end+=data_list[cnt_rank+1][se]
                    if cnt_rank+1<len(data_list):
                        start_new = sum(data_list[cnt_rank+1][:start])
                        end_new = start_new+sum_end
  
                    rank_next = cnt_rank+2
                    while rank_next<len(dis_shift_sim):
          
                        end_next = 0
         
                        for nse in range(start_new,end_new):
                    
                            if rank_next<len(data_list):
                                end_next+=data_list[rank_next][nse]
                            if dis_shift_sim[rank_next][nse] == rank_next - cnt_rank-1:
       
                                dis_bottom2root[rank_next][nse]+=rank_l
                        if rank_next<len(data_list):
                            start_next =sum(data_list[rank_next][:start_new])
                            start_new = start_next
                            end_new = start_new+end_next
                        rank_next+=1

                cnt_l+=1
                    
        cnt_rank+=1
    return dis_bottom2root

# Determine whether it is a subtree.
def is_subtree2(list1,list2):
    
    if len(list1)==1:

        if list1[0][0]<=list2[0][0]:
            return True,[0],[0]
        else:
            return False,[-1],[-1]
    else:# len != 1

        index_1 =list(range(len(list1[1]))) 
        
        cnt_head = 0  
        flag_find = False 
        list_head = [] # init None
        
        flag_not_find = False # true can not find subtree
        while len(list_head) < len(list1[1])  and flag_find==False and flag_not_find==False:
            index_2_used = [] # used node in list2 
            index_1_cnt = 0 # index of current node 
            flag_brak = False # break flag 
            for cnt_l1_1 in index_1: 
                if flag_brak:
                    break
                flag_12 = False # true: subtree, false: not subtree 
                l1_1 = list1[1][cnt_l1_1]
                s_list1 = spilt_list(list1,1,cnt_l1_1)
                
                cnt_l2_1 = 0
                
                while cnt_l2_1 < (len(list2[1])) and flag_12==False: # next layer of list2 
                    if cnt_l2_1 not in index_2_used:
                        l1_2 = list2[1][cnt_l2_1] 
                        if l1_2>= l1_1:
                            s_list2 = spilt_list(list2,1,cnt_l2_1)

                            if is_len(s_list1,s_list2):
                                

                                flag_12,_,_= is_subtree2(s_list1,s_list2)
                            else:
                                
                                flag_12 = False
                        else:
                            flag_12 = False
                        if flag_12==True: # find a subtree 
                            # this node used 
                            
                            index_2_used.append(cnt_l2_1)
 
                        
                    cnt_l2_1+=1
                if flag_12==False:
                    if cnt_l1_1 not in list_head and cnt_l1_1!=index_1[0]: 
                        
                        # reset index_1 
                        flag_brak = True # break
                        
                        
                        flag_2_bk = False
                        for used_2 in reversed(index_2_used):# find in used 
                            if flag_2_bk == True:
                                break
                            l1_2_used = list2[1][used_2]# cnt of used nodes 
                            if l1_2_used>= l1_1:
                                s_list2 = spilt_list(list2,1,used_2)
                                if is_len(s_list1,s_list2):
                                    flag_12_2,_,_= is_subtree2(s_list1,s_list2)
                                else:
                                     flag_12_2 = False
                            else:
                                 flag_12_2 = False
                            if flag_12_2==True:# break 
                                re_used2 = used_2 # used twice 
                                list_head.append(cnt_l1_1) 
   
                                index1_used2 = index_2_used.index(re_used2)
                                
                                l1_re = index_1[index1_used2]
             
                                index_1[index1_used2] = cnt_l1_1
                                
                                index_1[index_1_cnt] = l1_re
                                
                                flag_2_bk = True 
                            else:
                                return False,[-1],[-1]
                                    
                        
                    else: # this node can not be head 
                        return False,[-1],[-1]
                else:
                    index_1_cnt+=1
            if flag_12:
                flag_find = True
                
            cnt_head+=1
        
        return flag_find,index_1,index_2_used
# Find the optimal sublist.
def find_closest_index(lst, target, exclude_indices):
    
    closest_index = 0
    min_difference = 100  # init inf 
    cnt_i = 0
    for l in lst:
        if cnt_i not  in exclude_indices:
          # skip used index 
            difference = abs(l - target)  # calculate dis
            if difference < min_difference:
                min_difference = difference
                closest_index = cnt_i
        cnt_i+=1

    return closest_index,min_difference
# Find the optimal list dis
def get_best_list_index(list1,list2):# greedy
    
    ls_dis = 0
    
    flag_in_2 = True
    index_in_2 = [] 
    for l1 in list1:
        if l1 not in list2:
            flag_in_2 = False
            break
        else:
            index_in_2.append(list2.index(l1))
    if flag_in_2:
        
        return 0,index_in_2,flag_in_2
    else:
        
        exclude_indices = []
        for l1 in list1:
            
            a_index,a_dis = find_closest_index(list2, l1, exclude_indices)
           
            exclude_indices.append(a_index)
            ls_dis+=a_dis

        return float("{:.8f}".format(ls_dis)),exclude_indices,flag_in_2
# Identify the closest subtree.
def get_tree_dis_new_3(data_t1,dis_sp1,data_t2,dis_sp2,data2_name,cnt_it,shuffle=False,flag_dis=False,dis1_pre=[],dis_pre=[]):
    # dis_pre: pre dis, init [] 

    
    if len(data_t1)==1:# last layer
        
        if dis_sp2[0][0]>0:
            
            list_dis_add = dis_pre+dis_sp2[0]
        else:
            list_dis_add = copy.deepcopy(dis_pre)
        dis_2_com = [dsp+sum(list_dis_add) for dsp in dis_sp2[1]]
        if dis_sp1[0][0]>0:
            list_dis_add1 = dis1_pre+dis_sp1[0]
            
        else:
            list_dis_add1 = copy.deepcopy(dis1_pre)
        dis_1_com = [dsp+sum(list_dis_add1) for dsp in dis_sp1[1]]
 
        min_dis_2,choice_index,flag_in_12 = get_best_list_index(dis_1_com,dis_2_com)
        if flag_dis:
            
                list_cg_flag = [0,1]
                f_cg = random.choice(list_cg_flag)
                if f_cg==0:
                    # if len()
                    dis2_std = np.std(dis_2_com)
                    dis1_std = np.std(dis_1_com)
                    min_1_d = min(dis_1_com) - dis1_std
                    max_1_d = max(dis_1_com) + dis1_std
                    dis_2_new_minmax = [d2c for d2c in dis_2_com if d2c>=min_1_d and d2c<=max_1_d]
                    if len(dis_2_new_minmax)>=len(dis_1_com):
                        dis_1_com = list(random.sample(dis_2_new_minmax,len(dis_1_com)))
                    else:
                        for d1c_i,d1c_m in enumerate(dis_1_com):
                    
                            dis_change_rate_list = list(np.linspace(d1c_m-dis2_std,d1c_m+dis2_std,5))
                    
                            d1c_new = random.choice(dis_change_rate_list)
                            dis_1_com[d1c_i] = d1c_new
                   
                    min_dis_2,choice_index,_ = get_best_list_index(dis_1_com,dis_2_com)
            
            

        name_sp = [data2_name[1][ci] for ci in choice_index]
        name_list = [[data2_name[0][0]]]
        name_list.append(list(name_sp))
        
        return abs(float("{:.8f}".format(min_dis_2))),name_list 
    else:
       
        flag_subtree,index_t1,index_t2 = is_subtree2(data_t1,data_t2)
        
        if flag_subtree :
            
            if dis_sp2[0][0]>0:
                dis_pre = dis_pre+dis_sp2[0]
            if dis_sp1[0][0]>0:
                dis1_pre = dis1_pre+dis_sp1[0]
            
            min_dis_sum = 0 # sum of dis of multiple nodes
            name_list=[[data2_name[0][0]]] # upon node name 
            
            name_list_next = []
            flag_it2 = False 
            
            if flag_it2==False:
                
                min_dis_sum = 0 
                
                index_2_used = []
                dis_sp1_index = list(range(len(dis_sp1[1])))
                dis_sp2_index = list(range(len(dis_sp2[1])))
                if cnt_it==1:
                    if shuffle:
                        random.shuffle(dis_sp1_index)
                    
                        print("shuffle",dis_sp1_index)
                for dis1 in dis_sp1_index:
                    flag_find = False 
                    min_dis12 = 10000000
                    dis_sp1_split = split_dis_list(data_t1,1,dis1,dis_sp1)
                    
                    data_t1_split = spilt_list(data_t1,1,dis1)
                    
                    for dis2 in dis_sp2_index:
                        data_t2_split = spilt_list(data_t2,1,dis2)
                        
                        if dis2 not in index_2_used and is_subtree2(data_t1_split,data_t2_split)[0]:
                            dis_sp2_split = split_dis_list(data_t2,1,dis2,dis_sp2)
                            data_name2_split = split_dis_list(data_t2,1,dis2,data2_name)
                            
                            dis_120,name_list_return = get_tree_dis_new_3(data_t1_split,dis_sp1_split,
                            data_t2_split,dis_sp2_split,
                            data_name2_split,cnt_it+1,shuffle,flag_dis,dis1_pre.copy(),dis_pre.copy())
                            if dis_120<min_dis12:
                                flag_find = True 
                                min_dis12 = dis_120
                                index_2 = dis2
                                dis_2_best = name_list_return[:]
                                if min_dis12==0:
                                    break
                                
                                
                        
                    if flag_find: # find a subtree
                        index_2_used.append(index_2)
                        min_dis_sum+=min_dis12
                    
                        if len(name_list_next)>0:
                            name_list_next_new = [name_list_next[ni]+dis_2_best[ni] \
                                                  for ni in range(len(dis_2_best))]
                        
                        else:
                            name_list_next_new = dis_2_best
                        name_list_next = name_list_next_new.copy()
                        
                if flag_find:
                    name_list+=name_list_next
                    
                    return float("{:.8f}".format(min_dis_sum)),name_list
                else:
                    
                    min_dis_sum = 0
                
                    for id1,id2 in zip(index_t1,index_t2):
                        dis_sp1_split = split_dis_list(data_t1,1,id1,dis_sp1) 
                        
                        data_t1_split = spilt_list(data_t1,1,id1)
                        

                        data_t2_split = spilt_list(data_t2,1,id2)
                        dis_sp2_split = split_dis_list(data_t2,1,id2,dis_sp2)
                        data_name2_split = split_dis_list(data_t2,1,id2,data2_name)
                        
                        dis_12,name_list_return = get_tree_dis_new_3(data_t1_split,dis_sp1_split,
                                data_t2_split,dis_sp2_split,data_name2_split,cnt_it+1,                                                                     shuffle,flag_dis,dis1_pre.copy(),dis_pre.copy())
                        min_dis_sum+=dis_12
                        if len(name_list_next)>0:
                            name_list_next_new = [name_list_next[ni]+name_list_return[ni] \
                                                  for ni in range(len(name_list_return))]
                        else:
                            name_list_next_new = name_list_return[:]

                        name_list_next = name_list_next_new.copy()

                    name_list+=name_list_next
                    
                    return float("{:.8f}".format(min_dis_sum)),name_list

        else:
            print("no sbu tree")
            return 300000,["nosub"]
# Split the list. (name list)
def spilt_list(data,n,i):

    data_spilt = []
    num_n = data[n][i] 
    data_spilt.append([num_n])
    pre_n = sum(data[n][:i]) 
    
    cnt_n = 1# traverse next layers 
    while cnt_n+n<len(data):
        
        data_spilt.append(data[cnt_n+n][pre_n:pre_n+num_n])
        num_n_new = sum(data[cnt_n+n][pre_n:pre_n+num_n])
        num_n = num_n_new
        pre_n = sum(data[cnt_n+n][:pre_n]) 
        
        cnt_n+=1
    return data_spilt
# Split the dis list.
def split_dis_list(data,n,i,level):

    level_spilt = []
    num_n = data[n][i] 
    level_spilt.append([level[n][i]])
    pre_n = sum(data[n][:i]) 
    
    cnt_n = 1 
    while cnt_n+n<len(data):
        
        level_spilt.append(level[cnt_n+n][pre_n:pre_n+num_n])
        num_n_new = sum(data[cnt_n+n][pre_n:pre_n+num_n])
        num_n = num_n_new
        pre_n = sum(data[cnt_n+n][:pre_n]) 
        
        cnt_n+=1
    level_spilt.append(level[cnt_n+n][pre_n:pre_n+num_n]) 
    return level_spilt
# Find a sublist in List 2 where its elements are greater than or equal to those in List 1.
def is_len(list1,list2):
    for l1,l2 in zip(list1,list2):
        
        if len(l1) >len(l2):
            
            return False
        else:
            l1_sort = sorted(l1)            
            l2_sort = sorted(l2)
            index_l2_end = -1
            
            for i_1 in range(len(l1)):
                if l1_sort[index_l2_end]>l2_sort[index_l2_end]:
                    return False
                index_l2_end-=1
    return True

# List addition.
def my_add_list(list1,list2):
    list_all = []
    for l1,l2 in zip(list1,list2):
        list_all.append(list(l1+l2))
    return list_all

# sample a multiway tree (env mode)
def get_sp_cnt_data_level(path_meanshift,data_tree_train,dis_tree_train,
    shift_tree_train,position_phy,uni_name_phy_bac,cnt_name_phy_std_bac,uni_name_phy_ar,cnt_name_phy_std_ar):
    level = ["d","p","c","o","f","g"]
    list_generate = []
    l="d"
    path_x = os.path.join(path_meanshift,"x_"+l+".json")
    path_label = os.path.join(path_meanshift,"label_"+l+".json")
    x = np.loadtxt(path_x)
    label = np.loadtxt(path_label)
    # with open(path_x,"r") as f:
    #     x = json.load(f)
    # with open(path_label,"r") as f:
    #     label = json.load(f)

        
    d1 = int(random.choice((x[label==0])))
    d2 = int(random.choice((x[label==1])))
    

    domain=[d1,d2]
    
    domain_name = ["Archaea","Bacteria"]
    i_del = 0
    for d_v in domain:
        if d_v==0:
            del domain[i_del]
            del domain_name[i_del]
        else:
            i_del+=1
    list_generate.append(domain)

    cnt_p = 0
    
    if len(list_generate[0])>1:# two domians
        len_p11 = list_generate[0][0]
        len_p22 = list_generate[0][1]
        len_p1 = max(len_p11,len_p22)
        len_p2 = min(len_p11,len_p22)
        
        list_p1 = list(np.random.choice(uni_name_phy_bac,size =len_p1, p=cnt_name_phy_std_bac,replace=False))
        list_p2 = list(np.random.choice(uni_name_phy_ar,size =len_p2, p=cnt_name_phy_std_ar,replace=False))
        
        list_p = list_p1+list_p2
    else:
        len_p = sum(list_generate[0])
        list_p = np.random.choice(uni_name_phy_bac,size =len_p, p=cnt_name_phy_std_bac,replace=False)
    
    cnt_p_g = [] 
    p_choice_list = []
    for p in list_p:
        
        p_choice_all = random.choice(position_phy[p])
        p_choice = (p_choice_all[0],p_choice_all[1])
        p_choice_list.append(p_choice)
        data_p_choice = spilt_list(data_tree_train[p_choice[0]],1,p_choice[1])
        cnt_p_g.append(sum(data_p_choice[-1]))
    sorted_indices = sorted(range(len(cnt_p_g)), key=lambda i: cnt_p_g[i],reverse=True)
    
    for p_i in sorted_indices:
        p = list_p[p_i]
        
        p_choice = p_choice_list[p_i]

        data_p_choice = spilt_list(data_tree_train[p_choice[0]],1,p_choice[1])
        dis_p_choice = split_dis_list(data_tree_train[p_choice[0]],1,p_choice[1],dis_tree_train[p_choice[0]])
        shift_p_choice = split_dis_list(data_tree_train[p_choice[0]],1,p_choice[1],shift_tree_train[p_choice[0]])
        
        
        if cnt_p==0:
            
            data_p = data_p_choice
            
            dis_p = dis_p_choice
            shift_p = shift_p_choice
        else:
            data_p = my_add_list(data_p,data_p_choice)
            dis_p = my_add_list(dis_p,dis_p_choice)
            shift_p = my_add_list(shift_p,shift_p_choice)
            
        cnt_p+=1
    
    return [list_generate[0]]+data_p,[[0 for lg0 in  list_generate[0]]]+dis_p,[[0 for lg0 in  list_generate[0]]]+shift_p


# sample a multiway tree (sample mode)
def get_sp_cnt_data_level_sample(len_d_list,data_tree_train,dis_tree_train,
    shift_tree_train,position_phy,uni_name_phy_bac,cnt_name_phy_std_bac,uni_name_phy_ar,cnt_name_phy_std_ar):
    list_generate = []
    d1 = random.choice(len_d_list[0]) #ar
    d2 = random.choice(len_d_list[1]) #bac
    domain=[d1,d2]
    domain_name = ["Archaea","Bacteria"]
    i_del = 0
    for d_v in domain:
        if d_v==0:
            del domain[i_del]
            del domain_name[i_del]
        else:
            i_del+=1
    list_generate.append(domain)
    cnt_p = 0
    
    len_p1 = d1
    len_p2 = d2
    
    
    list_p1 = []
    list_p2 = []
    if len_p1>0:
        list_p1 = list(np.random.choice(uni_name_phy_ar,size =len_p1, p=cnt_name_phy_std_ar,replace=False))
    if len_p2>0:
        list_p2 = list(np.random.choice(uni_name_phy_bac,size =len_p2, p=cnt_name_phy_std_bac,replace=False))
    list_p = list_p1+list_p2
    
    cnt_p_g = [] 
    p_choice_list = []
    for p in list_p:
        p_choice_all = random.choice(position_phy[p])
        p_choice = (p_choice_all[0],p_choice_all[1])
        p_choice_list.append(p_choice)
        data_p_choice = spilt_list(data_tree_train[p_choice[0]],1,p_choice[1])
        cnt_p_g.append(sum(data_p_choice[-1]))
    sorted_indices = sorted(range(len(cnt_p_g)), key=lambda i: cnt_p_g[i],reverse=True)

    for p_i in sorted_indices:
        p = list_p[p_i]

        p_choice = p_choice_list[p_i]

        data_p_choice = spilt_list(data_tree_train[p_choice[0]],1,p_choice[1])
        dis_p_choice = split_dis_list(data_tree_train[p_choice[0]],1,p_choice[1],dis_tree_train[p_choice[0]])
        shift_p_choice = split_dis_list(data_tree_train[p_choice[0]],1,p_choice[1],shift_tree_train[p_choice[0]])

        if cnt_p==0:    
            data_p = data_p_choice
            
            dis_p = dis_p_choice
            shift_p = shift_p_choice
        else:
            data_p = my_add_list(data_p,data_p_choice)
            dis_p = my_add_list(dis_p,dis_p_choice)
            shift_p = my_add_list(shift_p,shift_p_choice)
            
        cnt_p+=1
    return [list_generate[0]]+data_p,[[0 for lg0 in  list_generate[0]]]+dis_p,[[0 for lg0 in  list_generate[0]]]+shift_p

