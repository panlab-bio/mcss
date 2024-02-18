# -*- coding: utf-8 -*-
# Simulate clustering on the statistically obtained data_all, merely learning from real data without actual simulation
path_meanshift = os.path.join(path_sp_cnt,"meanshift")
if not os.path.exists(path_meanshift):
    os.makedirs(path_meanshift)
# Save the results of mean shift.
def save_meanshift(data_all,path_meanshift = path_meanshift):
    level = ["d","p","c","o","f","g"]

    criteria = (cv2.TERM_CRITERIA_EPS+cv2.TERM_CRITERIA_MAX_ITER,10,1) 
    # Evaluation metrics for k-means 
    K = 2
    # k represents the number of lineages under this hierarchy.
    # The number of clusters obtained by mean shift is denoted as mk
    flag =True # Check if it's definitive; use k-means for boundaries.
    for index_l,i in enumerate(data_all):
        path_x = os.path.join(path_meanshift,"x_"+level[index_l]+".json")
        path_label = os.path.join(path_meanshift,"label_"+level[index_l]+".json")
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
        np.savetxt(path_x, x)
        np.savetxt(path_label, my_label)
        # with open(path_x,"w") as f:
        #     json.dump(x,f)
        # with open(path_label,"w") as f:
        #     json.dump(my_label,f)
save_meanshift(data_all)
