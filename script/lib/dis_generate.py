import numpy as np
import pandas as pd
import itertools
import subprocess
import os
import random
from scipy import stats
import pickle 
import copy

def get_distance_dict_node_alltree_2(n,merge_name,my_columns,dis_n):
    #返回层与层之间的dis，输入的是某一层
    clm=my_columns[n]
    len_my_columns = len(my_columns)
    dis_l = dict()#这个之前是存每个linege和层之间的dis,现在把他变成和node之间dis，key为lingeage，val为列表
    dis_list=[]
    #其实只需要判断一下lineage前面的node是不是实体节点?
    # dict_above = dict()#用于记录一个点是不是实体点
    list_above = [] #只要出现一个实体点,就把他放到该列表中,判断一个新加的点的上面的点是不是在其中
    cnt_n=1
    while (n+cnt_n)<len_my_columns:
        dis_list.append([])
        cnt_n+=1
    
    for i in dis_n[n].keys():#i是某个lineage,如细菌界

        list_dis_i = []#记录这个linege的node dis列表
        list_sp_i_node = [] #记录有上层实体点的sp
        cmd_lineage="merge_name[merge_name."+my_columns[n]+"=="+"i]"
        # 获得这个 lineage 的 dataframe
        lineage = eval(cmd_lineage)
        
        cnt_n = 1 #1时为下一层，2时为下下层
        #首先认为dis_n中的key都是实体点了,下面那个判断len(genus)其实有点多余
        # n+cnt_n从1开始，到哪结束？到dis_n[-1]结束，dis_n 0界 1门 2纲 3目4科5属而len_my_columns6
        if n==0 and len(lineage)==1:
            list_dis_i = [{},{},{},{},{}]
        else:
            while (n+cnt_n)<len_my_columns-1:
                list_dis_next = [] #linege i 与某一层node的dis
                list_dis_next_dict = dict()
                cmd_genus="lineage."+my_columns[n+cnt_n]+".unique()"

                genus = eval(cmd_genus)#genus是界下面的门,uni列表,后面要遍历所有的门,然后判断是不是实体点?
                # l_genus+=len(genus
                # if len(genus)>=1:#其实这个判断界是不是实体点,但是其实是多余的
                if len(genus)>1 or n==0:
                    #是个实体点
                    # if cnt_n>1: #不是临近层要判断是不是上面有实体点,eg界下属上面是不是有纲
                    #下下层 
                    for j in genus:#界下某个具体的门

                        if j in dis_n[n+cnt_n].keys() :#这个也是实体点
                            cmd_lineage_j="merge_name[merge_name."+my_columns[n+cnt_n]+"=="+"j]"
                            #这个是j这个具体的门对应的dataframe
                            lineage_j = eval(cmd_lineage_j)

                            cnt_below = cnt_n+1 #用于记录下一层的下下层
                            while (n+cnt_below)<=len_my_columns-1:

                                cmd_genus_j="lineage_j."+my_columns[n+cnt_below]+".unique()"
                                genus_j = eval(cmd_genus_j)#门下纲的uni列表
                                cnt_below+=1

                            #界下的门,这次没根据len判断门是不是实体点,因为默认在dis中的都是,如果他不是,那么到他作为根节点时也会报错
                            s1 = list(dis_n[n+cnt_n][j].keys())[0]
                            # print("n",n,"n+cnt_n",n+cnt_n,"超过了dis_n边界,cnt_n为",cnt_n)
                            a_dis = np.float('%.6f'%(dis_n[n][i][s1]-dis_n[n+cnt_n][j][s1]))

                            list_dis_next.append(a_dis)
                            list_dis_next_dict[j] = a_dis
                            dis_list[cnt_n-1].append(a_dis)

                list_dis_i.append(list_dis_next_dict)
                cnt_n+=1
        dis_i_sp = [ key for key,val in dis_n[n][i].items()]
        list_dis_next_sp = []
        list_dis_next_dict_sp = dict()
        for keysp in dis_i_sp:
            # if keysp not in list_below:
            a_dis_sp = dis_n[n][i][keysp]
            list_dis_next_sp.append(a_dis_sp)
            list_dis_next_dict_sp[keysp] = a_dis_sp
            dis_list[-1].append(a_dis_sp)
        list_dis_i.append(list_dis_next_dict_sp)
        dis_l[i] = list_dis_i
    return dis_l,dis_list

def get_distance2(i,merge_name,my_columns,taxonomy,path_nw,path_sim,abs_path_dir):
    dis = []
    #i是第i列,i从零开始
    cmd_unique = "merge_name."+my_columns[i]+".unique()"
    l_unique = eval(cmd_unique)#所有的门
    # print("Species unique个数",len(l_unique))
    cnt=0
    ss = 0
    for j in l_unique:#j是某个特定的门
        domain = merge_name[merge_name[my_columns[i]]==j][my_columns[0]].unique()[0]
        # print("domain",domain)
        cmd_select = "merge_name[merge_name."+my_columns[i] + "==j]"
        aa = eval(cmd_select)#exec不能返回值，用eval，某个门对应的dataframe        
        len_s = len(aa) #该门的物种数目
        cmd_len = "len(aa."+my_columns[i+1]+".unique())"
        len_aa = eval(cmd_len)#该门下纲数目
        ss+=len_s#所有门物种数总和
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
    
    # print("len>1",len(dis),"species",cnt,"all",ss)
    return dis
    

def get_subtree_new(a,domain,path_nw,path_sim,abs_path_dir):
    #返回的是字典，key为物种名，值为dis
    path_data = os.path.join(abs_path_dir,"data")
    path_tg = os.path.join(path_sim,"tmp_data","tmp_g.txt")
    path_tg_root = os.path.join(path_sim,"tmp_data")
    a.to_csv(path_tg,index=False,header=None) 
    #这个好像成功了,好像其实不用保存名字,还是要保存
    if domain=="Bacteria":
        tree_name = os.path.join(path_data,"gtdb_tree/bac_sp_new.tree")
    else:
        tree_name = os.path.join(path_data,"gtdb_tree/ar_sp_new.tree")
        
    #其实名字不该是_g，因为界门纲目都可以用该函数，不过无所谓了，就先用这个名字吧
    path_prune = os.path.join(path_nw,"nw_prune")
    cmd_sub = path_prune+" -v "+ tree_name+" -f tmp_g.txt >sub_g.txt"
    
    # print(cmd_sub)
    #不保存不行，还是老问题，a是个数组，多行，但是cmd只要了一行(第一行)
    subprocess.run(cmd_sub,shell=True,cwd=path_tg_root)
    #/bin/sh: 1: nw_prune: Permission denied 出现这种问题可以用绝对路径
    path_distance = os.path.join(path_nw,"nw_distance")
    cmd_dis = path_distance+" -n sub_g.txt"
    sp = subprocess.run(cmd_dis,shell=True,cwd=path_tg_root,
                        capture_output=True,encoding='utf-8')
    
    dis_g = sp.stdout.split()

    idx=0
    dis_dict=dict()
    for i in range(int(len(dis_g)/2)):
        # print(idx,dis_g[idx],dis_g[idx+1])
        dis_dict[dis_g[idx]]=np.float(dis_g[idx+1])
        idx+=2             
    return dis_dict
    

def get_distance_new(i,merge_name,my_columns,taxonomy,path_nw,path_sim,abs_path_dir):
    #返回字典，key是门，值为每个物种到该门的dis
    dis = dict()
    cmd_unique = "merge_name."+my_columns[i]+".unique()"
    l_unique = eval(cmd_unique)
    # print("Lineage with name unique个数",len(l_unique))
    cnt=0
    ss = 0
    for j in l_unique:
        domain = merge_name[merge_name[my_columns[i]]==j][my_columns[0]].unique()[0]
        # print("domain",domain)
        cmd_select = "merge_name[merge_name."+my_columns[i] + "==j]"
        aa = eval(cmd_select)#exec不能返回值，用eval
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
            if i==0 and len_aa==1:     #给界开个绿灯
                #其实是门的个数
                dis_g = dict()
                aa_s = aa.Species
                # print(my_columns[i],"可以在只有一个物种的情况下计算")
                cmd_select_tax = "taxonomy[taxonomy."+my_columns[i] + "==j]"
                aa_tax = eval(cmd_select_tax)
                aa_s_tax = aa_tax.Species
                dis_g_tax = get_subtree_new(aa_s_tax,domain,path_nw,path_sim,abs_path_dir)
                for aa_s_i in aa_s:
                    dis_g[aa_s_i] = dis_g_tax[aa_s_i]
                dis[str(j)] = dis_g
            
    # print("len>1",len(dis),"species",cnt,"all",ss)
    return dis


def get_num_new_2(merge_name):
    
    data = merge_name
    # data_null = data.isnull()
    # data = data.fillna(axis=1,method='ffill')
    # data[data_null] =data[data_null] +'_nan'
    clu = data.columns
    datad =[y for y in ([data[data[clu[i]]==x][clu[i+1]].nunique() \
                         for x in data[clu[i]].unique()] for i in range(len(clu)-1))]
    datad_name = [y for y in ([x \
                         for x in data[clu[i]].unique()] for i in range(len(clu)))]
    return datad,datad_name

def get_dis_list_alltree(df_merge,path_tax,path_nw,path_sim,abs_path_dir):
    #类似与sum_fun_dis
    #但是sum_fun_dis是计算层与层之间的dis，而这个是记录node和node之间的距离
#     sp_cut=pd.read_csv(path_name,header=None,usecols=[0],sep="\t")
#     sp_cut.columns=["Species"]
    
    taxonomy =pd.read_csv(path_tax,sep="\t")
    # df_merge = pd.merge(taxonomy,sp_cut,on=["Species"])
    my_columns = df_merge.columns
    # data_sp = get_num_new(taxonomy)
    data_sp = []
    dis_n=[]
    for i in range(6):
        # print("i",i)
        dis_i = get_distance_new(i,df_merge,my_columns,taxonomy,path_nw,path_sim,abs_path_dir)
        dis_n.append(dis_i)#sp和lieage距离
    # print("get_distance")
    dis_5 = get_distance2(5,df_merge,my_columns,taxonomy,path_nw,path_sim,abs_path_dir)
    
    dis_g_new = dis_n[5] #属要么没有实体点，有实体店，只能和计算距离，所以这个还是对的
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

def get_up_node_cnt(l1,n,taxonomy):
    clu = taxonomy.columns
    cnt_up = 1#先设置初始值，是几无所谓，只要>1能让循环跑通就行
    cnt_level = 0
    name_up = l1
    while n-1>=0 and cnt_up<=1:
        name_up = taxonomy[taxonomy[clu[n]] == name_up][clu[n-1]].unique()[0] #上一层的lineage
        # cnt_up = taxonomy[taxonomy[clu[n-1]] == name_up][clu[n]].nuinque()
        cnt_up = taxonomy[taxonomy[clu[n-1]] == name_up][clu[n]].nunique()
        n-=1
        cnt_level+=1
    # print(n)
    return cnt_level,name_up

def get_dis_tree(dis_all_dict_alltree,datad_name,data_sp_alltree,taxonomy):
    cnt_rank = 0
    dis_tree = []
    shift_tree = []
    # clu = taxonomy.columns
    clu = list(dis_all_dict_alltree.keys())
    for rank in datad_name:
        dis_tree_rank = []
        shift_tree_rank = []
        if cnt_rank==0:#界全赋值成0

            dis_tree_rank = [0 for r0 in rank]
            shift_tree_rank = [0 for r0 in rank]

        else:
            #用来衡量这个节点的最近上层实体点
            cnt_l = 0
            for l in rank:
                if cnt_rank<len(data_sp_alltree):#到属
                    if data_sp_alltree[cnt_rank][cnt_l]>1:
                        # print(data_sp_alltree[cnt_rank][cnt_l])
                        cnt_level,name_up = get_up_node_cnt(l,cnt_rank,taxonomy)
                        # print(clu[cnt_rank-cnt_level],name_up,l)
                        dis_l = dis_all_dict_alltree[clu[cnt_rank-cnt_level]][0][name_up][cnt_level-1][l]
                        # print(dis_l)
                        dis_tree_rank.append(dis_l)
                        shift_tree_rank.append(cnt_level-1) #cnt_level从1开始的，所以1代表上一层，shift_tree中0表示上一层，所以要-1

                    else:
                        cnt_level,name_up = get_up_node_cnt(l,cnt_rank,taxonomy)
                        dis_tree_rank.append(-1)
                        #对于非实体店，虽没有dis，但是有上层节点数，但是好像没有不用管？因为不会用到非实体点的shift?
                        shift_tree_rank.append(cnt_level-1)
                else:
                    cnt_level,name_up = get_up_node_cnt(l,cnt_rank,taxonomy)
                    # print(clu[cnt_rank-cnt_level],name_up,l,cnt_level-1)
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

def get_root_node(data_sim,dis_sp_sim,dis_shift_sim):
    # print(len(data_sim[0]))
    if len(data_sim[0])==1:#看看有几个界，
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
    # print(dis_list)
    for rank in dis_bottom2root[:-1]:#到属结束
        cnt_l = 0 #记录第几个门
        if cnt_rank>=1: #从门开始
            # print("rank",rank)
            for rank_l in rank:#遍历每个门
                
                if rank_l>0: #实体点

                    start = sum(data_list[cnt_rank][:cnt_l]) #cnt_1从0 开始 如果前面没有应该是看 :1的值？因为尾端娶不到1
                    
                    end = start+data_list[cnt_rank][cnt_l]
                    # print("start,end",start,end)
                    sum_end = 0
                    for se in range(start,end): #下一层
            
                        if dis_bottom2root[cnt_rank+1][se]>0:
           
                            dis_bottom2root[cnt_rank+1][se]+=rank_l #只是改了一层
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

def is_subtree2(list1,list2):
    #判断左边的列表是不是右边列表的子树
    if len(list1)==1:

        if list1[0][0]<=list2[0][0]:
            return True,[0],[0]
        else:
            return False,[-1],[-1]
    else:#长度不是1了

        index_1 =list(range(len(list1[1]))) 
        # index_2 = [] #把2也记住
        cnt_head = 0 #用于记录多少次
        flag_find = False #没找到一个合适的顺序，其实好像和flag12重复了
        list_head = [] #不用cnt_head了，cnt_head不对，因为012 -> 102 -> 如果在0卡住了，可能又回去了
        #认为初始值为空，这比初始值为0的好处是我们遍历完一次后可以把它放进去，而不是判断他的下一个
        flag_not_find = False #如果为true，说明确定找不到了，比如有好几个都抢占同一个节点，那你怎么分也分不过来
        while len(list_head) < len(list1[1])  and flag_find==False and flag_not_find==False:
            index_2_used = [] #记录list2中已经用过的点，每次换head都要重新清空
            index_1_cnt = 0 #要记录当前节点的位置，就是他的index虽然是cnt_l1_1，如2 ，他的位置可能是index_1_cnt 0
            flag_brak = False #这一组index_1有没有中途遇到问题
            for cnt_l1_1 in index_1: #[2],遍历list1下一层,1是下一层
                if flag_brak:
                    break
                flag_12 = False #默认是不是子树，返回true才认为是子树
                l1_1 = list1[1][cnt_l1_1]#[2]
                s_list1 = spilt_list(list1,1,cnt_l1_1)
                # print("list1 分割之后",s_list1)
                cnt_l2_1 = 0
                # print("现在找的是1中",cnt_l1_1)
                while cnt_l2_1 < (len(list2[1])) and flag_12==False: #遍历list2下一层
                    if cnt_l2_1 not in index_2_used:
                        l1_2 = list2[1][cnt_l2_1]#[2,2]中的2
                        if l1_2>= l1_1:
                            s_list2 = spilt_list(list2,1,cnt_l2_1)

                            if is_len(s_list1,s_list2):
                                # print("长度符合")
                                # print(s_list1,s_list2,"\n")

                                flag_12,_,_= is_subtree2(s_list1,s_list2)
                            else:
                                # print("长度不符合\n")
                                flag_12 = False
                        else:
                            flag_12 = False
                        if flag_12==True:#找到了一个子树
                            #这个点用过了，下面就不再找这个节点了
                            # print("当前used",index_2_used)
                            index_2_used.append(cnt_l2_1)
 
                        
                    cnt_l2_1+=1
                if flag_12==False:#
                    if cnt_l1_1 not in list_head and cnt_l1_1!=index_1[0]: 
                        
                        #重新设置index_1
                        flag_brak = True #中断了
                        # del index_1[index_1_cnt]
                        # index_1 = [cnt_l1_1]+index_1
                        #现在是换两个冲突的
                        flag_2_bk = False
                        for used_2 in reversed(index_2_used):#从这些用过的里面找
                            if flag_2_bk == True:
                                break
                            l1_2_used = list2[1][used_2]# 用过的节点的数目信息
                            if l1_2_used>= l1_1:
                                s_list2 = spilt_list(list2,1,used_2)
                                if is_len(s_list1,s_list2):
                                    flag_12_2,_,_= is_subtree2(s_list1,s_list2)
                                else:
                                     flag_12_2 = False
                            else:
                                 flag_12_2 = False
                            if flag_12_2==True:#别忘了break
                                re_used2 = used_2 #这个被使用了两次，这是2的节点，要看上次他被使用时的1节点是啥
                                list_head.append(cnt_l1_1) #这个节点换过一次了
   
                                index1_used2 = index_2_used.index(re_used2)
                                #1节点下标对应的1节点
                                l1_re = index_1[index1_used2]
             
                                index_1[index1_used2] = cnt_l1_1
                                #新的换成旧的
                                index_1[index_1_cnt] = l1_re
                                # print("更换后的1是",index_1)
                                flag_2_bk = True #要更换了，下面的不用看了
                            else:
                                return False,[-1],[-1]
                                    
                        
                    else:#这个节点当过一次头，还找不到最优解，那没办法了
                        return False,[-1],[-1]
                else:
                    index_1_cnt+=1
            if flag_12:
                flag_find = True
                
            cnt_head+=1
        # print(flag_find,index_1)
        return flag_find,index_1,index_2_used

def find_closest_index(lst, target, exclude_indices):
    # print("exclude_indices",exclude_indices)
    closest_index = 0
    min_difference = 100  # 初始差值设置为正无穷大
    cnt_i = 0
    for l in lst:
        if cnt_i not  in exclude_indices:
          # 跳过指定的下标
            difference = abs(l - target)  # 计算当前元素与目标值的差值
            if difference < min_difference:
                min_difference = difference
                closest_index = cnt_i
        cnt_i+=1

    return closest_index,min_difference

def get_best_list_index(list1,list2):#贪心
    # print("--")
    ls_dis = 0
    # exclude_indices = [-1]
    flag_in_2 = True
    index_in_2 = [] 
    for l1 in list1:
        if l1 not in list2:
            flag_in_2 = False
            break
        else:
            index_in_2.append(list2.index(l1))
    if flag_in_2:
        # print("完全找到")
        return 0,index_in_2,flag_in_2
    else:
        # print("不能完全找到")
        exclude_indices = []
        for l1 in list1:
            # print(l1)
            a_index,a_dis = find_closest_index(list2, l1, exclude_indices)
            # print(a_index,a_dis)
            exclude_indices.append(a_index)
            ls_dis+=a_dis

        return float("{:.8f}".format(ls_dis)),exclude_indices,flag_in_2
def get_tree_dis_new_3(data_t1,dis_sp1,data_t2,dis_sp2,data2_name,cnt_it,shuffle=False,flag_dis=False,dis1_pre=[],dis_pre=[]):
    #dis_pre用来记录前面的dis 列表，最开始的时候为[] 
    #输入的是树的距离列表，应该是level列表，因为要判断是不是0
    #level是只到属水平的，但是dis是到种水平的，所以要改一下level列表
    # print("dis列表",dis_sp1,"---",dis_sp2)
    # print("data列表",data_t1,"---",data_t2)
    # print("\n")
    
    if len(data_t1)==1:#到最后一层节点了
        # print()
        # min_dis_1,choice_index = get_best_list_index(dis_sp1[-1],dis_sp2[-1])#没改dis之前的距离
        # shift_pre+=dis_shift_2[0]
        # shift1_pre+=dis_shift_1[0]
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
            # for d1c_i,d1c_m in enumerate(dis_1_com):
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
                    # # print(d1c_m,dis2_std,d1c_m-dis2_std,d1c_m+dis2_std)
                    # # dis_change_rate_list = [0.0002,-0.0002,0,0,0,0.0005,-0.0005,0.001,-0.001]
                            dis_change_rate_list = list(np.linspace(d1c_m-dis2_std,d1c_m+dis2_std,5))
                    # dis_change_rate_list = list(np.linspace(d1c_m-dis1_std,d1c_m+dis1_std,100))
                    # # print(dis_change_rate_list)
                            d1c_new = random.choice(dis_change_rate_list)
                            dis_1_com[d1c_i] = d1c_new
                    # print(d1c_m,"d1c_new",d1c_new)
                    # print()        
                # else:
                #     print("d1c_m",d1c_m)
                #     print()
                
                # print(dis_change_rate_list)
                
                # d2c_list = []
                # d1c_m = np.mean(dis_1_com)
                # for d2c_i,d2c_v in enumerate(dis_2_com):
                #     if abs(d2c_v - d1c_m) <= d1c_m*0.01:
                #         d2c_list.append(d2c_i)
                # if len(d2c_list)>= len(dis_1_com):
                #     choice_index = random.sample(d2c_list,len(dis_1_com))
                    
                # dis_change_rate = random.choice(dis_change_rate_list)
                # dis_1_com[d1c_index] = dis_1_com[d1c_index]+dis_1_com[d1c_index]*dis_change_rate
                    min_dis_2,choice_index,_ = get_best_list_index(dis_1_com,dis_2_com)
            
            

        name_sp = [data2_name[1][ci] for ci in choice_index]
        name_list = [[data2_name[0][0]]]
        name_list.append(list(name_sp))
        # print("最后一层的头节点的名字",name_list)
            
        # print("min_dis_2,choice_index",min_dis_2,choice_index,"\n")
        # print(min_dis_1,min_dis_2)
        return abs(float("{:.8f}".format(min_dis_2))),name_list #min_dis_2是正确的dis
    else:
        # flag_subtree = False#判断是不是sub
        # print(data_t1,data_t2)
        flag_subtree,index_t1,index_t2 = is_subtree2(data_t1,data_t2)
        # if flag_subtree:
        #     print("data_t1是data_t2的子树")
        #     # print(data_t1,data_t2)
        # else:
        #     print("data_t1不是data_t2的子树")
            # print(data_t1,data_t2)

        # cnt_sub1 = 0 #记录dis_t1第几个子树
        if flag_subtree :#到时候看一下flag_subtree的时候返回值是啥
            # if cnt_it==0:
            #     shift_pre=[]
            #     dis_pre=[]
            # shift_pre+=dis_shift_2[0]
            # shift1_pre+=dis_shift_1[0]
            if dis_sp2[0][0]>0:
                dis_pre = dis_pre+dis_sp2[0]
            if dis_sp1[0][0]>0:
                dis1_pre = dis1_pre+dis_sp1[0]
            # print("改变之后的dis_pre",dis_pre)
            min_dis_sum = 0 #记录下面多个点的dis和
            name_list=[[data2_name[0][0]]] #先把最上面的点的名字加上
            # print("头节点名字",name_list)
            name_list_next = []
            flag_it2 = False #默认情况下不走index_t2，但是第一条走不通了再走
            # if cnt_it<=10:
            
                
                # print()
                # print(list_my_p)
            if flag_it2==False:
                # if cnt_it<2:
                #     print("选项1")
                min_dis_sum = 0 
                # print("cnt_it",cnt_it)
                # print(data_t1,data_t2)
                # print(data_t1)
                # cnt_sub1 = 0
                index_2_used = []
                dis_sp1_index = list(range(len(dis_sp1[1])))
                dis_sp2_index = list(range(len(dis_sp2[1])))
                if cnt_it==1:
                    if shuffle:
                        random.shuffle(dis_sp1_index)
                    
                        print("shuffle",dis_sp1_index)
                for dis1 in dis_sp1_index:
                    flag_find = False #一开始这个1没有找到合适的2
                    min_dis12 = 10000000
                    dis_sp1_split = split_dis_list(data_t1,1,dis1,dis_sp1)
                    # dis_shift_1_split = split_dis_list(data_t1,1,cnt_sub1,dis_shift_1) 
                    #因为dis_shift_1和dis_sp1形状一样
                    data_t1_split = spilt_list(data_t1,1,dis1)
                    # name_1_split = split_dis_list(data_t1,1,cnt_sub1,name_1) 
                    # cnt_sub2 = 0
                    for dis2 in dis_sp2_index:
                        data_t2_split = spilt_list(data_t2,1,dis2)
                        # print(data_t1_split,"\n",data_t2_split,"\n")
                        if dis2 not in index_2_used and is_subtree2(data_t1_split,data_t2_split)[0]:
                            dis_sp2_split = split_dis_list(data_t2,1,dis2,dis_sp2)
                            data_name2_split = split_dis_list(data_t2,1,dis2,data2_name)
                            # dis_shift_2_split = split_dis_list(data_t2,1,cnt_sub2,dis_shift_2)
                            dis_120,name_list_return = get_tree_dis_new_3(data_t1_split,dis_sp1_split,
                            data_t2_split,dis_sp2_split,
                            data_name2_split,cnt_it+1,shuffle,flag_dis,dis1_pre.copy(),dis_pre.copy())
                            if dis_120<min_dis12:
                                flag_find = True #其实只要进循环了，就认为找到了，不过这个更好一点
                                min_dis12 = dis_120
                                index_2 = dis2
                                dis_2_best = name_list_return[:]
                                if min_dis12==0:
                                    break
                                # print("name_list_return",name_list_return)
                                
                        # cnt_sub2+=1
                    if flag_find: #找到了再赋值
                        index_2_used.append(index_2)
                        min_dis_sum+=min_dis12
                    # print("index_2_used",index_2_used,"min_dis_sum",min_dis_sum)
                        if len(name_list_next)>0:
                            name_list_next_new = [name_list_next[ni]+dis_2_best[ni] \
                                                  for ni in range(len(dis_2_best))]
                        # print("name_list_next_new",name_list_next_new)
                        else:
                            name_list_next_new = dis_2_best
                        name_list_next = name_list_next_new.copy()
                        # cnt_sub1+=1
                if flag_find:
                    name_list+=name_list_next
                    # if cnt_it<2:
                    #     print(cnt_it,data2_name[0][0],"找到了")
                    #     print(name_list[:3])
                    
                # if min_dis_sum>10000:
                #     print("小cnt",cnt_it)
                    # print("贪心",min_dis_sum,name_list)
                    return float("{:.8f}".format(min_dis_sum)),name_list
                else:
                    # if cnt_it<2:
                    #     print("选项2")
                    # if len(set(name_1[-1])&set(list_my_p[1]))>0 and len(set(datad_name)&set(list_my_p[1]))>0:
                    #     print("222")
                    # if cnt_it<2:
                    #     print(cnt_it,data2_name[0][0],"这次没找到，换")
                    min_dis_sum = 0
                # print("cnt_it",cnt_it)
                # print("index_t1,index_t2",index_t1,index_t2)
                # print(data_t1,data_t2,"\n")
                    for id1,id2 in zip(index_t1,index_t2):
                        dis_sp1_split = split_dis_list(data_t1,1,id1,dis_sp1) 
                        #1 cnt_sub1，下一层的第cnt_sub1个节点的子树
                        # dis_shift_1_split = split_dis_list(data_t1,1,id1,dis_shift_1) 
                        #因为dis_shift_1和dis_sp1形状一样
                        data_t1_split = spilt_list(data_t1,1,id1)
                        # name_1_split = split_dis_list(data_t1,1,id1,name_1) 

                        data_t2_split = spilt_list(data_t2,1,id2)
                        dis_sp2_split = split_dis_list(data_t2,1,id2,dis_sp2)
                        data_name2_split = split_dis_list(data_t2,1,id2,data2_name)
                        # dis_shift_2_split = split_dis_list(data_t2,1,id2,dis_shift_2)
                        # print("data_t1_split,data_t2_split",data_t1_split,data_t2_split)
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
                    # print("锚定",min_dis_sum,name_list)
                    return float("{:.8f}".format(min_dis_sum)),name_list

        else:
            print("no sbu tree")
            return 300000,["nosub"]

def spilt_list(data,n,i):
#把[[2],[1,1],[2,1],[2,2,2]]拆分，第n层第i个
    data_spilt = []
    num_n = data[n][i] #他下面有几个
    data_spilt.append([num_n])
    pre_n = sum(data[n][:i])#他前面有几个
    # print("下一层前面有",pre_n,"这一层下面有",num_n,"列表",data_spilt)
    cnt_n = 1#用于遍历下面的层
    while cnt_n+n<len(data):
        # print("获得的是",data[cnt_n+n],"第",pre_n,"到",pre_n,"+",num_n,"个",data[cnt_n+n][pre_n:pre_n+num_n])
        data_spilt.append(data[cnt_n+n][pre_n:pre_n+num_n])
        num_n_new = sum(data[cnt_n+n][pre_n:pre_n+num_n])
        num_n = num_n_new
        pre_n = sum(data[cnt_n+n][:pre_n])#pre是下一层的，所以要后面求
        # print("num_n",num_n,"pre_n",pre_n,pre_n+num_n,data[cnt_n+n],data[cnt_n+n][pre_n:pre_n+num_n])
        # print("下一层前面有",pre_n,"这一层下面有",num_n,"列表",data_spilt)
        cnt_n+=1
    return data_spilt

def split_dis_list(data,n,i,level):
#把[[2],[1,1],[2,1],[2,2,2]]拆分，第n层第i个
    level_spilt = []
    num_n = data[n][i] #他下面有几个
    level_spilt.append([level[n][i]])
    pre_n = sum(data[n][:i])#他前面有几个
    # print("下一层前面有",pre_n,"这一层下面有",num_n,"列表",data_spilt)
    cnt_n = 1#用于遍历下面的层
    while cnt_n+n<len(data):
        # print("获得的是",data[cnt_n+n],"第",pre_n,"到",pre_n,"+",num_n,"个",data[cnt_n+n][pre_n:pre_n+num_n])
        level_spilt.append(level[cnt_n+n][pre_n:pre_n+num_n])
        num_n_new = sum(data[cnt_n+n][pre_n:pre_n+num_n])
        num_n = num_n_new
        pre_n = sum(data[cnt_n+n][:pre_n])#pre是下一层的，所以要后面求
        # print("num_n",num_n,"pre_n",pre_n,pre_n+num_n,data[cnt_n+n],data[cnt_n+n][pre_n:pre_n+num_n])
        # print("下一层前面有",pre_n,"这一层下面有",num_n,"列表",data_spilt)
        cnt_n+=1
    level_spilt.append(level[cnt_n+n][pre_n:pre_n+num_n]) #比level多了最后一层
    return level_spilt

def is_len(list1,list2):
    for l1,l2 in zip(list1,list2):
        # print(l1,l2)
        if len(l1) >len(l2):
            # print("l1比l2长了",l1,l2)
            return False
        else:
            l1_sort = sorted(l1)            
            l2_sort = sorted(l2)
            index_l2_end = -1
            # while l2_sort[index_l2_strat] < l1_sort[0]:
            #     index_l2_strat+=1
            for i_1 in range(len(l1)):
                if l1_sort[index_l2_end]>l2_sort[index_l2_end]:
                    return False
                index_l2_end-=1
    return True
def my_add_list(list1,list2):
    list_all = []
    for l1,l2 in zip(list1,list2):
        list_all.append(list(l1+l2))
    return list_all

def get_sp_cnt_data_level(path_meanshift,data_tree_train,dis_tree_train,
    shift_tree_train,position_phy,uni_name_phy_bac,cnt_name_phy_std_bac,uni_name_phy_ar,cnt_name_phy_std_ar):
    level = ["d","p","c","o","f","g"]
    list_generate = []
    l="d"
    path_x = os.path.join(path_meanshift,"x_"+l+".pkl")
    path_label = os.path.join(path_meanshift,"label_"+l+".pkl")

    with open(path_x,"rb") as f:
        x = pickle.load(f)
    with open(path_label,"rb") as f:
        label = pickle.load(f)

        # print(x[label==0])
    # print(x[label==1])
    d1 = random.choice(x[label==0])
    d2 = random.choice(x[label==1])
    # if min(d1,d2)>1:
    #     print(x[label==0])
    #     print(x[label==1])

    domain=[d1,d2]
    # print("domain",domain)
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
    # len_p = len(list_generate[-1])
    if len(list_generate[0])>1:#有古菌
        len_p11 = list_generate[0][0]
        len_p22 = list_generate[0][1]
        len_p1 = max(len_p11,len_p22)
        len_p2 = min(len_p11,len_p22)
        
        list_p1 = list(np.random.choice(uni_name_phy_bac,size =len_p1, p=cnt_name_phy_std_bac,replace=False))#选出来的具体门
        list_p2 = list(np.random.choice(uni_name_phy_ar,size =len_p2, p=cnt_name_phy_std_ar,replace=False))#选出来的具体门
        # print(list_p1,list_p2)
        list_p = list_p1+list_p2
    else:
        len_p = sum(list_generate[0])
        list_p = np.random.choice(uni_name_phy_bac,size =len_p, p=cnt_name_phy_std_bac,replace=False)
    
    cnt_p_g = [] #记录门下属的个数
    p_choice_list = []
    for p in list_p:#先选出门，根据选的门的个数来生成
        #或者反正我就几个门，能不能门的最优解？
        #先把大门放前面
        p_choice_all = random.choice(position_phy[p])
        p_choice = (p_choice_all[0],p_choice_all[1])
        p_choice_list.append(p_choice)
        data_p_choice = spilt_list(data_tree_train[p_choice[0]],1,p_choice[1])
        cnt_p_g.append(sum(data_p_choice[-1]))
    sorted_indices = sorted(range(len(cnt_p_g)), key=lambda i: cnt_p_g[i],reverse=True)
    # print("选择的p",list_p,p_choice_list)
    # print("sorted_indices",sorted_indices)
    for p_i in sorted_indices:
        p = list_p[p_i]
        # p_choice_all = random.choice(position_phy[p])

        # p_choice = (p_choice_all[0],p_choice_all[1])
        p_choice = p_choice_list[p_i]

        data_p_choice = spilt_list(data_tree_train[p_choice[0]],1,p_choice[1])
        dis_p_choice = split_dis_list(data_tree_train[p_choice[0]],1,p_choice[1],dis_tree_train[p_choice[0]])
        shift_p_choice = split_dis_list(data_tree_train[p_choice[0]],1,p_choice[1],shift_tree_train[p_choice[0]])
        
        # print("data_p_choice",data_p_choice)
        if cnt_p==0:
            
            data_p = data_p_choice
            # print(data_p)
            dis_p = dis_p_choice
            shift_p = shift_p_choice
        else:
            data_p = my_add_list(data_p,data_p_choice)
            dis_p = my_add_list(dis_p,dis_p_choice)
            shift_p = my_add_list(shift_p,shift_p_choice)
            # print(data_p)
        cnt_p+=1
    # print(list_generate[0])
    return [list_generate[0]]+data_p,[[0 for lg0 in  list_generate[0]]]+dis_p,[[0 for lg0 in  list_generate[0]]]+shift_p



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
    # print("len_p1,len_p2",len_p1,len_p2)
    
    list_p1 = []
    list_p2 = []
    if len_p1>0:
        list_p1 = list(np.random.choice(uni_name_phy_ar,size =len_p1, p=cnt_name_phy_std_ar,replace=False))#选出来的具体门
    if len_p2>0:
        list_p2 = list(np.random.choice(uni_name_phy_bac,size =len_p2, p=cnt_name_phy_std_bac,replace=False))#选出来的具体门
    list_p = list_p1+list_p2
    
    cnt_p_g = [] #记录门下属的个数
    p_choice_list = []
    for p in list_p:#先选出门，根据选的门的个数来生成
        #或者反正我就几个门，能不能门的最优解？
        #先把大门放前面
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
            # print(data_p)
            dis_p = dis_p_choice
            shift_p = shift_p_choice
        else:
            data_p = my_add_list(data_p,data_p_choice)
            dis_p = my_add_list(dis_p,dis_p_choice)
            shift_p = my_add_list(shift_p,shift_p_choice)
            # print(data_p)
        cnt_p+=1
    return [list_generate[0]]+data_p,[[0 for lg0 in  list_generate[0]]]+dis_p,[[0 for lg0 in  list_generate[0]]]+shift_p

