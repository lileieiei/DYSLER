
import numpy as np
import pandas as pd
import time
import os 

from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.preprocessing import StandardScaler

from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel
from sklearn.feature_selection import mutual_info_regression,mutual_info_classif
from sklearn.ensemble import GradientBoostingClassifier as GBC
from sklearn.model_selection import cross_validate, KFold
from sklearn.feature_selection import SequentialFeatureSelector
from sklearn.metrics import accuracy_score as ACC 
from sklearn.metrics import log_loss as logloss 
import matplotlib.pyplot as plt


from telcoFunc import *

import features_creation as fc
from features_creation import *
from tqdm import tqdm
import gc
from sklearn.feature_selection import RFE 
from sklearn.ensemble import GradientBoostingClassifier as GBC 
from sklearn.model_selection import cross_validate, KFold
import hyperopt
from hyperopt import hp, fmin, tpe, Trials, partial
from hyperopt.early_stop import no_progress_loss

data = pd.read_csv("./6_1data_exp_shifted.csv",index_col=0)
data_test=pd.read_csv('./6_1data_exp_test_shifted.csv',index_col=0)
features = data.drop(columns='group').copy()
labels = data['group'].copy()
scaler = StandardScaler()
data_keep=scaler.fit_transform(features)
ff=pd.DataFrame(data_keep,columns=features.columns)
MIC = mutual_info_classif(ff, labels, random_state=1412)
MIC = pd.Series(MIC, index=ff.columns)
MIC =MIC.sort_values(ascending=False)
MIC=MIC.index[:-6]
ff_mic=ff.loc[:,MIC.values] 

score_keep=[]
for i in range(10):
    print(i)
    param_grid_simple = {'n_estimators': hp.quniform("n_estimators",30,100,1)
                  ,"lr": hp.quniform("learning_rate",0.3,0.8,0.01)
                  ,"criterion": hp.choice("criterion",['friedman_mse', 'squared_error'])
                  ,"loss":hp.choice("loss",[ 'log_loss',"exponential"])
                  ,"max_depth": hp.quniform("max_depth",2,7,1)
                  ,"subsample": hp.quniform("subsample",0.2,0.9,0.01)
                  ,"max_features": hp.choice("max_features",["sqrt","auto","log2",2,3,4,5,7,8,9])
                  ,"min_impurity_decrease":hp.quniform("min_impurity_decrease",0,2,0.01)
                 }
    def hyperopt_objective(params):
        reg = GBC(n_estimators = int(params["n_estimators"])
                  ,learning_rate = params["lr"]
                  ,criterion = params["criterion"]
                  ,loss = params["loss"]
                  ,max_depth = int(params["max_depth"])
                  ,max_features = params["max_features"]
                  ,subsample = params["subsample"]
                  ,min_impurity_decrease = params["min_impurity_decrease"]
                  ,random_state=1412
                  ,verbose=False)
        validation_loss = cross_validate(reg,X_train_temp,labels
                                         ,scoring="roc_auc"
                                         ,cv=5
                                         ,verbose=False
                                         ,n_jobs=-1
                                         ,error_score='raise')
        return -np.mean(abs(validation_loss["test_score"]))
    def param_hyperopt(max_evals=100):

        #保存迭代过程
        trials = Trials()

        #设置提前停止
        early_stop_fn = no_progress_loss(400)

        #定义代理模型
        params_best = fmin(hyperopt_objective
                           , space = param_grid_simple
                           , algo = tpe.suggest
                           , max_evals = max_evals
                           , verbose=False
                           , trials = trials
                           , early_stop_fn = early_stop_fn
                          )

        #打印最优参数，fmin会自动打印最佳分数
        #print("\n","\n","best params: ", params_best,
         #     "\n")
        return params_best, trials
    # 创建容器
    rfe_res_search1 = []
    rfe_rs1_cv = []
    # 执行循环
    for i in tqdm(range(37)):

        i = 37 - i


        # 首次循环时，创建X_train_temp
        if i == 37:
            X_train_temp = (ff_mic).copy()    

        # 训练模型，然后带入RFE评估器
        params_best, trials = param_hyperopt(1500)
        mol=GBC(n_estimators=int(params_best['n_estimators'])
                ,learning_rate=params_best['learning_rate']
                ,loss=["deviance","exponential"][params_best['loss']]
                ,max_depth=int(params_best['max_depth'])
                ,max_features=["sqrt","auto","log2",1,2,3,4,5,7,8,9][params_best['max_features']]
                ,min_impurity_decrease=params_best['min_impurity_decrease']
                ,subsample=params_best['subsample']
                ,criterion=['friedman_mse', 'squared_error'][params_best['criterion']]
                ,random_state=1412)
        rfe_search = RFE(estimator=mol, n_features_to_select=i).fit(X_train_temp, labels)
        score = cross_validate(mol,X_train_temp,labels
                                         ,scoring="roc_auc"
                                         ,cv=5
                                         ,verbose=False
                                         ,n_jobs=-1
                                         ,error_score='raise') 
        X_train_temp = ff_mic[rfe_search.get_feature_names_out()]    

        # 搜索本轮被淘汰的特征，并记入rfe_res_search1
        print(np.mean(abs(score["test_score"])))
        rfe_rs1_cv.append(np.mean(abs(score["test_score"])))
        rfe_res_search1.append(rfe_search.feature_names_in_[rfe_search.ranking_ != 1])

    # 清除临时变量
    score_keep.append(rfe_rs1_cv[::-1])
    gc.collect()
    
data=pd.DataFrame(score_keep)
dd=data.T
dd.to_csv('rfe.csv')
