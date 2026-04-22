#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 13:47:19 2025

@author: swainasish
"""

#%% import libs
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial import distance_matrix
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import Ridge,RidgeCV
from sklearn.decomposition import MiniBatchDictionaryLearning,PCA,TruncatedSVD
import time
from scipy.stats import pearsonr,spearmanr
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RepeatedKFold
#%% supporting functions
def normalisation(df):
    df=pd.DataFrame(df)
    co = df.columns
    std = StandardScaler(with_mean=False)
    df = std.fit_transform(df)
    return pd.DataFrame(df,columns=co)
def weighted_concate_mat(df,x_cor,y_cor):
    cor_mat = np.array([x_cor,y_cor],dtype=np.float32).T
    tree = cKDTree(cor_mat)    
    distances, nei_indexes = tree.query(cor_mat, k=5)
    n_spot = df.shape[0]
    neighbour_mean_exp = pd.DataFrame([df.iloc[nei_indexes[i,:],:].mean().values for i in range(n_spot)],
                                      columns =df.columns )
    # for i in range(3):
    #     norm_val = 1/(i+2)
    #     print(norm_val)
    #     new_nei = pd.DataFrame([neighbour_mean_exp.iloc[nei_indexes[i,:],:].mean().values for i in range(n_spot)],
    #                                       columns =neighbour_mean_exp.columns )
    #     neighbour_mean_exp = (1-norm_val) * neighbour_mean_exp + norm_val*new_nei

    return neighbour_mean_exp

dict_learn = TruncatedSVD(n_components=50)
# dict_learn = PCA(n_components=50)
def small_sampeling(scdf,spdf,n_cell=500):
    n_sc,n_sp = scdf.shape[0],spdf.shape[0]
    rand_sc = np.random.choice(n_sc,n_cell)
    rand_sp = np.random.choice(n_sp,n_cell)
    scdf_subset = scdf.iloc[rand_sc,:]
    spdf_subset = spdf.iloc[rand_sp,:]
    sc_sp_cooncate = pd.concat([scdf_subset,spdf_subset])

    return sc_sp_cooncate

def reference_selection(sp_ann_df,sc_ann_df,class_list,sc_gt):
    class_list=np.array(class_list)
    class_names = np.unique(class_list)
    classwise_dfs = [sc_ann_df.iloc[np.where(class_list==i)[0],:] for i in class_names]
    #finding common genes 
    cg = sp_ann_df.columns
    for df in classwise_dfs:
        cg = np.intersect1d(cg, df.columns)
    sp_ann_df = sp_ann_df.loc[:,cg]
    classwise_dfs = [i.loc[:,cg] for i in classwise_dfs]
    #finding gene coexprsession patterns 
    target_corr = sp_ann_df.corr()
    reference_corr = [classwise_dfs[i].corr() for i in range(len(classwise_dfs))]
    #diag = 0
    target_corr[target_corr==1]=0
    target_corr[target_corr<0]=0
    for i in range(len(reference_corr)):
        tempdf = reference_corr[i]
        tempdf[tempdf==1]=0
        tempdf[tempdf<0]=0
        tempdf[np.isnan(tempdf)]=0
        reference_corr[i] = tempdf
    #calculate pcc 
    pearson_corr_dic={}
    pearson_corr_dic_median={}
    for i in range(len(reference_corr)):
        dt_name = class_names[i]
        tempdf = reference_corr[i]
        pr = np.array([pearsonr(target_corr.iloc[:,i],tempdf.iloc[:,i]).statistic for i in range(target_corr.shape[1])])
        pr[np.isnan(pr)]=0
        pearson_corr_dic[dt_name]=pr
        pearson_corr_dic_median[dt_name]=np.median(pr)
    pearson_corr_dic = pd.DataFrame(pearson_corr_dic)
    sns.violinplot(pearson_corr_dic)
    plt.title("Pearson's corr of spatial dataset to the scRNA-seq references")
    plt.show()
    #choosing dataset
    max_corr = np.max(list(pearson_corr_dic_median.values()))
    limit_val = max_corr-(max_corr*0.25)
    best_indexes =  np.where(list(pearson_corr_dic_median.values())>limit_val)[0]
    grp_names = np.array(list(pearson_corr_dic_median.keys()))[best_indexes]
    select_cells = [True if i in grp_names else False for i in class_list]
    
    final_dtset = sc_ann_df[select_cells]
    sc_gt = sc_gt[select_cells]
    final_dtset=final_dtset.reset_index(drop=True)
    sc_gt=sc_gt.reset_index(drop=True)
    print(f"SPBOOST_OUT: Dataset '{grp_names}' choosen based on PCC to the spatial dataset")
    return final_dtset,sc_gt
    
#%% 
def GeneEnhancement(scdf,spdf,genes_to_predict,x_cor,y_cor,reference_selection=False,batch_info =False):
    t1= time.time()
    scdf=scdf.reset_index(drop=True)
    spdf=spdf.reset_index(drop=True)
    common_genes = np.intersect1d(scdf.columns, spdf.columns)
    print(f"{len(common_genes)} number of common genes found.")
    sc_g = scdf.loc[:,genes_to_predict]
    # sc_g = normalisation(sc_g)
    # cord_mat = np.array([x_cor,y_cor],dtype=np.float32).T

    #slicing both df
    scdf_cg=scdf.loc[:,common_genes]
    spdf_cg=spdf.loc[:,common_genes]
    if reference_selection ==True:
        scdf_cg,sc_g = reference_selection(spdf_cg,scdf_cg,batch_info,sc_g)
    
    scdf_cg = normalisation(scdf_cg)
    spdf_cg = normalisation(spdf_cg)

    
    n_sc=scdf.shape[0]
    n_sp=spdf.shape[0]
    n_all=n_sc+n_sp
    
    # #calculate neighbour expression 
    
    sp_neighbour = weighted_concate_mat(spdf_cg , x_cor, y_cor)
    
   
    knn_reg = KNeighborsRegressor(n_neighbors=10)
    knn_reg.fit(spdf_cg,sp_neighbour)
    sc_neighbour = pd.DataFrame(knn_reg.predict(scdf_cg),columns=scdf_cg.columns)
    
    #concate-weighted-mat 
    sc_wei = pd.concat([scdf_cg,sc_neighbour],axis=1)
    sp_wei = pd.concat([spdf_cg,sp_neighbour],axis=1)
    
    #dict_learnig
    n_sample = np.min([n_sc,n_sp])
    combine_sc_sp = small_sampeling(sc_wei,sp_wei,n_cell=n_sample)
    dict_learn.fit(combine_sc_sp)
    print("DICT FIT DONE")
    sc_wei_dict = dict_learn.transform(sc_wei)
    sp_wei_dict = dict_learn.transform(sp_wei)
    # sc_wei_dict = sc_wei
    # sp_wei_dict = sp_wei
    

    #elastic net 
    # cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=1)
    # ratios = arange(0, 1, 0.01)
    # alphas = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.0, 1.0, 10.0, 100.0]
    # model = ElasticNetCV(l1_ratio=ratios, alphas=alphas, cv=cv, n_jobs=-1)
    
    # #predict the expression 
    ridge_cv = RidgeCV(alphas=[1,10,100,1000,10000,100000])
    ridge_cv.fit(sc_wei_dict,sc_g)
    #print(ridge_cv.alpha_)
    model = Ridge(alpha=ridge_cv.alpha_)
    model.fit(sc_wei_dict,sc_g)
    sp_pred = model.predict(sp_wei_dict)
    sp_pred[sp_pred<0]=0
    
    
    sp_pred = pd.DataFrame(sp_pred,columns=sc_g.columns)
    # sp_pred = normalisation(sp_pred)
    
    t2=time.time()
    print(f"Took {t2-t1}  Seconds.")
    return  sp_pred