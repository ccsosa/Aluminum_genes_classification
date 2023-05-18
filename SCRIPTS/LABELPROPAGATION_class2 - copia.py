
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
#import seaborn as sns
# In[3]:
filename = "D:/TESIS_PHD/RICENETPPI/BIG_COMP_W.emb"
#filename = "/users/ccsosaa/pecanpy/BIG_COMP_W.emb"
data = pd.read_csv(filename,header=None, delim_whitespace=True,index_col=0)

filename = "D:/TESIS_PHD/RICENETPPI/JJG_DEG/CANDIDATES_1.csv"
#filename = "/users/ccsosaa/pecanpy/CANDIDATES_COMBINED.csv"
candidates = pd.read_csv(filename,index_col=0,header=None)

X = data
y = candidates.iloc[:,:]

idx = y.index
# #SCALING

scaler = preprocessing.StandardScaler()
X2 = X.copy()
scaler.fit(X2)
X2 = pd.DataFrame(scaler.transform(X2))

import numpy as np
X2["X1"] = np.array(y.iloc[:,0])
X2.index = idx
X3 = X2[X2["X1"]!=0]
X4 = X3
X3_to = X3.loc[:, X2.columns!='X1']
y3_to = X3.loc[:, X2.columns=='X1']

X4 = X2[X2["X1"]==0]
y4 = X4["X1"]
X4 = X4.loc[:, X4.columns!='X1']

X_train, X_test, y_train, y_test = train_test_split(X3_to, y3_to, test_size=0.10, random_state = 0)

from sklearn.semi_supervised import LabelPropagation
label_prop_model = LabelPropagation(n_jobs=-1, max_iter=10000)
LP = label_prop_model.fit(X_train, y_train)

y_pred_LP = label_prop_model.predict(X_test)
label_prop_model.score(X_test, y_test)

from sklearn import metrics
fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred_LP, pos_label=1)
print(metrics.auc(fpr, tpr))


#Calculo de las metricas Precision, Recall, F1-score 
from sklearn.metrics import classification_report
report_LP = classification_report(y_test,y_pred_LP)
print(report_LP)

#training
#y_pred_LP2 = label_prop_model.predict(X_train)

#from sklearn.metrics import classification_report
#report_LP = classification_report(y_train,y_pred_LP2)
#print(report_LP)

prob_LP = label_prop_model.predict_proba(X4)
pred_LP = label_prop_model.predict(X4)
#print(prob_svm)

results = y4.copy()
results=pd.DataFrame(results)

results[["prob_1",]] = prob_LP
results[["pred_1"]] = prob_LP
results["final"] =0

for i in range(len(results["prob_1"])):
    if results["pred_1"][i] >= 0.9:
        results["pred_1"][i] = 1
    else:
        results["pred_1"][i] = 0


for i in range(len(results["final"])):
    if results["pred_1"][i] ==1:
        results["final"][i] = 1
    else:
        results["final"][i] = 0
     
print(len(results["final"][results["final"]==0]))
print(len(results["final"][results["final"]==1]))


#results.to_csv("D:/TESIS_PHD/RICENETPPI/out_pred.csv")
results.to_csv("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/RESULTS/WEIGHTED/out_pred_class1LP_JJG.csv")


# In[4]:
filename = "D:/TESIS_PHD/RICENETPPI/BIG_COMP.emb"
#filename = "/users/ccsosaa/pecanpy/BIG_COMP_W.emb"
data = pd.read_csv(filename,header=None, delim_whitespace=True,index_col=0)

filename = "D:/TESIS_PHD/RICENETPPI/JJG_DEG/CANDIDATES_1.csv"
#filename = "/users/ccsosaa/pecanpy/CANDIDATES_COMBINED.csv"
candidates = pd.read_csv(filename,index_col=0,header=None)

X = data
y = candidates.iloc[:,:]

idx = y.index
# #SCALING

scaler = preprocessing.StandardScaler()
X2 = X.copy()
scaler.fit(X2)
X2 = pd.DataFrame(scaler.transform(X2))

import numpy as np
X2["X1"] = np.array(y.iloc[:,0])
X2.index = idx
X3 = X2[X2["X1"]!=0]
X4 = X3
X3_to = X3.loc[:, X2.columns!='X1']
y3_to = X3.loc[:, X2.columns=='X1']

X4 = X2[X2["X1"]==0]
y4 = X4["X1"]
X4 = X4.loc[:, X4.columns!='X1']

X_train, X_test, y_train, y_test = train_test_split(X3_to, y3_to, test_size=0.10, random_state = 0)

from sklearn.semi_supervised import LabelPropagation
label_prop_model = LabelPropagation(n_jobs=-1, max_iter=10000)
LP = label_prop_model.fit(X_train, y_train)

y_pred_LP = label_prop_model.predict(X_test)
label_prop_model.score(X_test, y_test)

from sklearn import metrics
fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred_LP, pos_label=1)
print(metrics.auc(fpr, tpr))


#Calculo de las metricas Precision, Recall, F1-score 
from sklearn.metrics import classification_report
report_LP = classification_report(y_test,y_pred_LP)
print(report_LP)

#training
#y_pred_LP2 = label_prop_model.predict(X_train)

#from sklearn.metrics import classification_report
#report_LP = classification_report(y_train,y_pred_LP2)
#print(report_LP)

prob_LP = label_prop_model.predict_proba(X4)
pred_LP = label_prop_model.predict(X4)
#print(prob_svm)

results = y4.copy()
results=pd.DataFrame(results)

results[["prob_1",]] = prob_LP
results[["pred_1"]] = prob_LP
results["final"] =0

for i in range(len(results["prob_1"])):
    if results["pred_1"][i] >= 0.9:
        results["pred_1"][i] = 1
    else:
        results["pred_1"][i] = 0


for i in range(len(results["final"])):
    if results["pred_1"][i] ==1:
        results["final"][i] = 1
    else:
        results["final"][i] = 0
     
print(len(results["final"][results["final"]==0]))
print(len(results["final"][results["final"]==1]))


#results.to_csv("D:/TESIS_PHD/RICENETPPI/out_pred.csv")
results.to_csv("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/RESULTS/UNWEIGHTED/out_pred_class1LP_JJGs.csv")

