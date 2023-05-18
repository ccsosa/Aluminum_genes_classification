#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
#from sklearn.model_selection import train_test_split
#from sklearn.linear_model import LogisticRegressionCV
#from sklearn.metrics import accuracy_score
#from gensim.models import Word2Vec
from sklearn import preprocessing
import sklearn.metrics as sk_metrics
#import seaborn as sns

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


#filename = "D:/TESIS_PHD/RICENETPPI/BIG_COMP_W.emb"
filename = "/users/ccsosaa/pecanpy/BIG_COMP_W.emb"
data = pd.read_csv(filename,header=None, delim_whitespace=True,index_col=0)

#filename = "D:/TESIS_PHD/RICENETPPI/CANDIDATES_1.csv"
filename = "/users/ccsosaa/pecanpy/CANDIDATES_1.csv"
candidates = pd.read_csv(filename,index_col=0,header=None)

from pulearn import BaggingPuClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression


# In[76]:


#X = data.iloc[1:400,:]
#y = candidates.iloc[1:400,:]

X = data
y = candidates.iloc[:,:]

#print(X)


# In[78]:


#print(y)


# #SCALING

# In[8]:


#pu_estimator.fit(X, y)
#scaler = preprocessing.MinMaxScaler()
scaler = preprocessing.StandardScaler()
X2 = X.copy()
scaler.fit(X2)
X2 = pd.DataFrame(scaler.transform(X2))

# In[10]:

print("SVM")

from pulearn import BaggingPuClassifier
from sklearn.svm import SVC
estimator = SVC(C=10, kernel='rbf', gamma=0.4, probability=True)


# In[11]:


pu_estimator1 = BaggingPuClassifier(
    base_estimator=estimator,
    n_estimators = 10000,  # 10000 trees as usual
    n_jobs = -1,           # Use all cores
    max_samples = sum(y.iloc[:,0]), # Each training sample will be balanced
    bootstrap=True
)

# In[13]:


pu_estimator1.fit(X2, y.iloc[:,0])


# In[14]:


y_pred_svm = pu_estimator1.predict(X2)
pu_estimator1.score(X2, y.iloc[:,0])

from sklearn import metrics
fpr, tpr, thresholds = metrics.roc_curve(y, y_pred_svm, pos_label=1)
print(metrics.auc(fpr, tpr))

# In[17]:


#Calculo de las metricas Precision, Recall, F1-score 
from sklearn.metrics import classification_report
report_svm = classification_report(y,y_pred_svm)
print(report_svm)


# In[19]:


prob_svm = pu_estimator1.predict_proba(X2)[:,1]
#print(prob_svm)


# In[22]:

fixed_svm_m = fixed_metrics(y,y_pred_svm)
print("AUC","AUCcr","ACCURACY","BALANCED ACCURACY","F","MCC","ALPHA")
print(fixed_svm_m)


# In[23]:

# #LDA

# In[25]:

print("LDA")

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from pulearn import BaggingPuClassifier
#estimator = SVC(C=10, kernel='rbf', gamma=0.4, probability=True)
#estimator = LogisticRegression(random_state=0,max_iter=500)
estimator = LinearDiscriminantAnalysis()

pu_estimator_LDA = BaggingPuClassifier(
    base_estimator=estimator,
    n_estimators = 1000,  # 10000 trees as usual
    #max_samples = sum(y), # Balance the positives and unlabeled in each bag
    max_samples = sum(y.iloc[:,0]), # Each training sample will be balanced
    n_jobs = -1,           # Use all cores
    bootstrap=True
)


# In[26]:


pu_estimator_LDA.fit(X2, y.iloc[:,0])


# In[27]:


y_pred_LDA = pu_estimator_LDA.predict(X2)
pu_estimator_LDA.score(X2, y.iloc[:,0])


# In[28]:


from sklearn import metrics
fpr, tpr, thresholds = metrics.roc_curve(y, y_pred_LDA, pos_label=1)
print(metrics.auc(fpr, tpr))


# In[29]:


#Calculo de las metricas Precision, Recall, F1-score 
from sklearn.metrics import classification_report
report_lda = classification_report(y,y_pred_LDA)
print(report_lda)


# In[30]:


prob_lda = pu_estimator_LDA.predict_proba(X2)[:,1]


fixed_LDA_m = fixed_metrics(y,y_pred_LDA)
print("AUC","AUCcr","ACCURACY","BALANCED ACCURACY","F","MCC","ALPHA")
print(fixed_LDA_m)

# In[31]:


print("RANDOM FOREST")

from pulearn import BaggingPuClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
pu_estimator_RF = BaggingPuClassifier(
    base_estimator= RandomForestClassifier(), 
    #RandomForestClassifier(n_estimators=1000, bootstrap= 1000),#
    n_estimators = 10000,  # 10000 trees as usual
    max_samples = sum(y.iloc[:,0]), # Each training sample will be balanced
    #max_samples = sum(y), # Balance the positives and unlabeled in each bag
    n_jobs = -1,           # Use all cores
    bootstrap=True
)


# In[32]:


pu_estimator_RF.fit(X2, y.iloc[:,0])


# In[33]:


y_pred_RF = pu_estimator_RF.predict(X2)
pu_estimator_RF.score(X2, y.iloc[:,0])


# In[34]:


from sklearn import metrics
fpr, tpr, thresholds = metrics.roc_curve(y, y_pred_RF, pos_label=1)
print(metrics.auc(fpr, tpr))


# In[35]:


prob_rf = pu_estimator_RF.predict_proba(X2)[:,1]
#print(prob_svm)


# In[37]:


#Calculo de las metricas Precision, Recall, F1-score 
from sklearn.metrics import classification_report
report_rf = classification_report(y,y_pred_RF)
print(report_rf)


# In[74]:


#prob_rf
fixed_RF_m = fixed_metrics(y,y_pred_RF)
print("AUC","AUCcr","ACCURACY","BALANCED ACCURACY","F","MCC","ALPHA")
print(fixed_RF_m)

# In[39]:


results = y.copy()
results=pd.DataFrame(results)


# In[40]:


y_pred_svm


# In[70]:


results["prob_svm"] = prob_svm
results["prob_lda"] = prob_lda
results["prob_rf"] = prob_rf
results["pred_svm"] = prob_svm
results["pred_lda"] = prob_lda
results["pred_rf"] = prob_rf
results["pred_A_svm"] = y_pred_svm
results["pred_A_lda"] = y_pred_LDA
results["pred_A_rf"] = y_pred_RF


# In[71]:

#print(results)
fixed_met_df = pd.DataFrame(fixed_svm_m)# columns=["AUC","ACCURACY","BALANCED ACCURACY","F","MCC","ALPHA"])
fixed_met_df = fixed_met_df.append(fixed_RF_m, ignore_index=False)
fixed_met_df = fixed_met_df.append(fixed_LDA_m, ignore_index=True)
fixed_met_df.columns =["AUC","AUCcr","ACCURACYcr","BALANCED ACCURACYcr","Fcr","MCCcr","ALPHA"]
fixed_met_df['Algorithm'] = ['SVM', 'RF', 'LDA']

#fixed_met_df.to_csv("D:/TESIS_PHD/RICENETPPI/fixed_metrics.csv")
fixed_met_df.to_csv("/users/ccsosaa/pecanpy/fixed_metrics.csv")

# In[72]:


for i in range(len(results["pred_svm"])):
    if results["pred_svm"][i] >= 0.9:
        results["pred_svm"][i] = 1
    else:
        results["pred_svm"][i] = 0

for i in range(len(results["pred_lda"])):
    if results["pred_lda"][i] >= 0.9:
        results["pred_lda"][i] = 1
    else:
        results["pred_lda"][i] = 0

for i in range(len(results["pred_rf"])):
    if results["pred_rf"][i] >= 0.9:
        results["pred_rf"][i] = 1
    else:
        results["pred_rf"][i] = 0


# In[77]:


#print(results)


# In[65]:

#results.to_csv("D:/TESIS_PHD/RICENETPPI/out_pred.csv")
results.to_csv("/users/ccsosaa/pecanpy/out_pred.csv")

