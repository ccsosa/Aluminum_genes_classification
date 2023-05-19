# -*- coding: utf-8 -*-
"""
Created on Fri May 19 11:39:46 2023

@author: cami_
"""


#!/usr/bin/env python
# coding: utf-8

# In[1]:
#import matplotlib.pyplot as plt
#import numpy as np
import pandas as pd
#import networkx as nx
#import csv
#from sys import stdin
#from networkx.algorithms.community import greedy_modularity_communities
#from networkx.algorithms.community import k_clique_communities
#from networkx.algorithms.community import asyn_lpa_communities
#from networkx.algorithms.community import modularity
#from cdlib import algorithms, readwrite, viz
#from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split
#from sklearn.linear_model import LogisticRegressionCV
#from sklearn.metrics import accuracy_score
#from gensim.models import Word2Vec
from sklearn import preprocessing
import sklearn.metrics as sk_metrics

from sklearn.semi_supervised import LabelPropagation

from pulearn import BaggingPuClassifier
#import seaborn as sns

# In[2]:
#################
def fixed_metrics(y,y_pred):
    from sklearn import metrics
    import pandas as pd
    
    

    cm = sk_metrics.confusion_matrix(y, y_pred)

    results_partial = y.copy()
    results_partial=pd.DataFrame(results_partial)
    results_partial["pred"] = y_pred
    results_partial = results_partial[(results_partial.iloc[:,0]<1)]
    
    alfa = sum(results_partial["pred"])/results_partial.shape[0]
    beta = 1
    
    fpr, tpr, thresholds = metrics.roc_curve(y, y_pred, pos_label=1)

    
    TN = cm[0,0]
    FN = cm[1,0]
    TP = cm[1,1]
    FP = cm[0,1]
    
    theta = (TP+FP)/(TP+FN+TN+FP)
    n=FP/(TN+FP)
    #pi = (TP+FN)/(TP+FN+TN+FP)
    c= len(y[(y.iloc[:,0]==1)])/(len(y[(y.iloc[:,0]==1)])+len(y[(y.iloc[:,0]==0)]))
    gamma = TP/(TP+FN)
    if beta-alfa==0:
        gamma_cr = float('inf')
        n_cr = float('inf')
    else:
        gamma_cr = ((beta-alfa)**-1)*((1-alfa)*gamma-(1-beta)*n)
        n_cr = ((beta-alfa)**-1)*((beta*n)-(alfa*gamma))
        
    pi_cr = (c*beta)+((1-c)*alfa)
    ##metrics
    ACC_cr = (pi_cr*gamma_cr) + ((1-pi_cr)*(1-n_cr))
    bacc_cr = (1 + (gamma_cr - n_cr))/2
    f_cr = (2*pi_cr*gamma_cr)/(pi_cr+theta)
    mcc_cr = (pi_cr*(1-pi_cr)/theta*(1-theta))*(gamma_cr-n_cr)
    
    auc_r = metrics.auc(fpr, tpr)
    auc_cr = (auc_r-(1-(beta-alfa)))/(beta-alfa)
    alfa = sum(results_partial["pred"])/results_partial.shape[0]
    
    metrics_df =[[auc_r,auc_cr,ACC_cr,bacc_cr,f_cr,mcc_cr,alfa]]
    return(metrics_df)

#################
    
# In[3]:
filename = "/users/ccsosaa/pecanpy//BIG_COMP_W.emb"
#filename = "/users/ccsosaa/pecanpy/BIG_COMP_W.emb"
data = pd.read_csv(filename,header=None, delim_whitespace=True,index_col=0)

filename = "/users/ccsosaa/pecanpy/CANDIDATES_1.csv"
#filename = "/users/ccsosaa/pecanpy/CANDIDATES_COMBINED.csv"
candidates = pd.read_csv(filename,index_col=0,header=None)

X = data
y = candidates.iloc[:,:]
#y[y==0] = -1
idx = y.index

scaler = preprocessing.StandardScaler()
X2 = X.copy()
scaler.fit(X2)
X2 = pd.DataFrame(scaler.transform(X2))


label_prop_model = LabelPropagation(n_jobs=-1, max_iter=10000)

pu_estimator1 = BaggingPuClassifier(
    base_estimator=label_prop_model,
    n_estimators = 10000,  # 10000 trees as usual
    n_jobs = -1,           # Use all cores
    max_samples = sum(y.iloc[:,0]), # Each training sample will be balanced
    bootstrap=True
)

pu_estimator1.fit(X2, y.iloc[:,0])

y_pred_LP = pu_estimator1.predict(X2)
pu_estimator1.score(X2, y.iloc[:,0])

from sklearn import metrics
fpr, tpr, thresholds = metrics.roc_curve(y, y_pred_LP, pos_label=1)
print(metrics.auc(fpr, tpr))


prob_LP = pu_estimator1.predict_proba(X2)[:,1]

fixed_LP_m = fixed_metrics(y,y_pred_LP)
print("AUC","AUCcr","ACCURACY","BALANCED ACCURACY","F","MCC","ALPHA")
print(fixed_LP_m)


results = y.copy()
results=pd.DataFrame(results)
results["prob_LP"] = prob_LP
results["pred_LP"] = prob_LP
results["pred_A_LP"] = y_pred_LP


for i in range(len(results["pred_LP"])):
    if results["prob_LP"][i] >= 0.9:
        results["pred_LP"][i] = 1
    else:
        results["pred_LP"][i] = 0

results.to_csv("/users/ccsosaa/pecanpy/weighted/out_pred_LAB_PROP.csv")


# In[4]:
    
filename = "/users/ccsosaa/pecanpy/BIG_COMP.emb"
#filename = "/users/ccsosaa/pecanpy/BIG_COMP_W.emb"
data = pd.read_csv(filename,header=None, delim_whitespace=True,index_col=0)

filename = "/users/ccsosaa/pecanpy/CANDIDATES_1.csv"
#filename = "/users/ccsosaa/pecanpy/CANDIDATES_COMBINED.csv"
candidates = pd.read_csv(filename,index_col=0,header=None)

X = data
y = candidates.iloc[:,:]
#y[y==0] = -1
idx = y.index

scaler = preprocessing.StandardScaler()
X2 = X.copy()
scaler.fit(X2)
X2 = pd.DataFrame(scaler.transform(X2))


label_prop_model = LabelPropagation(n_jobs=-1, max_iter=10000)

pu_estimator1 = BaggingPuClassifier(
    base_estimator=label_prop_model,
    n_estimators = 10000,  # 10000 trees as usual
    n_jobs = -1,           # Use all cores
    max_samples = sum(y.iloc[:,0]), # Each training sample will be balanced
    bootstrap=True
)

pu_estimator1.fit(X2, y.iloc[:,0])

y_pred_LP = pu_estimator1.predict(X2)
pu_estimator1.score(X2, y.iloc[:,0])

from sklearn import metrics
fpr, tpr, thresholds = metrics.roc_curve(y, y_pred_LP, pos_label=1)
print(metrics.auc(fpr, tpr))


prob_LP = pu_estimator1.predict_proba(X2)[:,1]

fixed_LP_m = fixed_metrics(y,y_pred_LP)
print("AUC","AUCcr","ACCURACY","BALANCED ACCURACY","F","MCC","ALPHA")
print(fixed_LP_m)


results = y.copy()
results=pd.DataFrame(results)
results["prob_LP"] = prob_LP
results["pred_LP"] = prob_LP
results["pred_A_LP"] = y_pred_LP


for i in range(len(results["pred_LP"])):
    if results["prob_LP"][i] >= 0.9:
        results["pred_LP"][i] = 1
    else:
        results["pred_LP"][i] = 0

results.to_csv("/users/ccsosaa/pecanpy/unweighted/out_pred_LAB_PROP.csv")