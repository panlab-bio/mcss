# -*- coding: utf-8 -*-
#对统计得到的data_all进行聚类模拟，只是学习真实数据，但是没模拟
path_meanshift = os.path.join(path_sp_cnt,"meanshift")
if not os.path.exists(path_meanshift):
    os.makedirs(path_meanshift)
#将meanshif的结果保存起来
def save_meanshift(data_all,path_meanshift = path_meanshift):
    level = ["d","p","c","o","f","g"]

    criteria = (cv2.TERM_CRITERIA_EPS+cv2.TERM_CRITERIA_MAX_ITER,10,1) 
    #kmeans的评价指标
    K = 2
    #k代表了这个层级下面lineage的个数
    #那么我meanshif得到的聚类数目是mk （meanshift K)
    flag =True #看看是不是界，界用kmeans
    for index_l,i in enumerate(data_all):
        path_x = os.path.join(path_meanshift,"x_"+level[index_l]+".pkl")
        path_label = os.path.join(path_meanshift,"label_"+level[index_l]+".pkl")
        if flag:
            x = np.array(i).reshape(-1,1)
            Z = np.array(i)
            Z = np.float32(Z)
            ret,my_label,center = cv2.kmeans(Z,K,criteria,None,10,cv2.KMEANS_PP_CENTERS)
            flag=False
        else:
            x = np.array(i).reshape(-1,1)
            n_samples=min(2000,int(len(i)))
            
            bandwidth = estimate_bandwidth(x, quantile=0.3,n_jobs=32,n_samples=n_samples)
            if bandwidth == 0:
                bandwidth = 0.2
            clustering = MeanShift(bandwidth=bandwidth, bin_seeding=True,n_jobs=32).fit(x)
            my_label = clustering.labels_
        print(path_x,path_label)
        with open(path_x,"wb") as f:
            pickle.dump(x,f)
        with open(path_label,"wb") as f:
            pickle.dump(my_label,f)
save_meanshift(data_all)