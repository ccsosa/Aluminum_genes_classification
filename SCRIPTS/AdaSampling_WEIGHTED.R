require(caret);require(data.table)
require(leaps);require(Metrics);require(AdaSampling)
###############################################################################s
################################################################################
#https://www.analyticsvidhya.com/blog/2021/04/smote-and-best-subset-selection-for-linear-regression-in-r/
dir <- "/users/ccsosaa/pecanpy"
################################################################################
#loading x and y
#x <- as.data.frame(fread(paste0(dir,"/","BIG_COMP_s.emb"),header = F))#,sep = "\t")
x <- as.data.frame(fread(paste0(dir,"/weighted/","BIG_COMP_W.emb"),header = F))#,sep = "\t")
row.names(x) <- x$V1; x$V1 <- NULL
y <- as.data.frame(fread(paste0(dir,"/","CANDIDATES_1.csv"),header = F))#,sep = "\t")
#y <- as.data.frame(fread(paste0(dir,"/JJG_DEG/","CANDIDATES_1.csv"),header = F))#,sep = "\t")
row.names(y) <- y$V1; y$V1 <- NULL    

y[,1] <- factor(y[,1],levels = c("1","0"))
y$names <- row.names(y)

################################################################################
#x <- x[1:500,]

#y <- y[y$names %in% row.names(x),]

set.seed(1234)
################################################################################
#https://towardsdatascience.com/building-classifiers-with-biased-classes-adasampling-comes-to-the-rescue-8212814264e3
#splitting in positive and negative
Ps <- y$names[which(y$V2=="1")]
Ns <- y$names[which(y$V2=="0")]
#y$names <- NA
################################################################################
#SCALING IT
x_scaled <- scale(x,center = T,scale = T)

################################################################################
# Apply AdaSampling method on the noisy data
#using all data
preds <- AdaSampling::adaSample(Ps, Ns, train.mat=x_scaled, test.mat=x_scaled, classifier = "lda",s=1,C = 10000)
preds_2 <- AdaSampling::adaSample(Ps, Ns, train.mat=x_scaled, test.mat=x_scaled, classifier = "svm",s=1,C = 10000)
preds_3 <- AdaSampling::adaSample(Ps, Ns, train.mat=x_scaled, test.mat=x_scaled, classifier = "logit",s=1,C = 10000)
#preds_4 <- AdaSampling::adaSample(Ps, Ns, train.mat=x_scaled, test.mat=x_scaled, classifier = "knn",s=1,C = 1000)

################################################################################
preds_list <- list(preds,preds_2,preds_3)#,,preds_4)
# accuracies <- list()
preds_list2 <- list() 
# preds_ensembl <- data.frame(matrix(ncol=7,nrow=nrow(preds_list[[1]])))
# preds_ensembl[,1] <- row.names(preds_list[[1]])


cutoff <- 0.9
for(i in 1:length(preds_list)){
 #i <- 1
  
  preds_1_fix <- data.frame(gene=row.names(preds_list[[i]]),prob=preds_list[[i]][,1],class=NA)
  for(j in 1:nrow(preds_1_fix)){

    if(preds_list[[i]][j,1]>=cutoff & preds_list[[i]][j,2]<cutoff){
      preds_1_fix$class[[j]] <- "0"
    } else if(preds_list[[i]][j,2]>=cutoff & preds_list[[i]][j,1]<cutoff){
      preds_1_fix$class[[j]] <- "1"
    } else if(preds_list[[i]][j,1]>=cutoff & preds_list[[i]][j,2]>=cutoff) {
      preds_1_fix$class[[j]] <- NA
    } else {
      preds_1_fix$class[[j]] <- "0"
    }
    
  };rm(j)
  preds_1_fix$class <- factor(preds_1_fix$class,levels = c("1","0"))
  if(i == 1){
    colnames(preds_1_fix) <- c("genes","prob_lda","pred_lda")
  } else   if(i == 2){
    colnames(preds_1_fix) <- c("genes","prob_svm","pred_svm")
  }  else if(i == 3){
    colnames(preds_1_fix) <- c("genes","prob_logit","pred_logit")
  }
  
  preds_list2[[i]] <- preds_1_fix
  # accuracies[[i]] <- caret::confusionMatrix(data = preds_1_fix[,3], reference = y[,1],mode = "everything", positive="1")
  # preds_ensembl[,i+1] <- as.numeric(as.character(preds_1_fix[,3]))
  
  };rm(i)


# preds_ensembl[,6] <- rowSums(preds_ensembl[,2:5],na.rm = T)
# 
# for(i in 1:nrow(preds_ensembl)){
#   if(preds_ensembl[i,6]>=2){
#     preds_ensembl[i,6] <- "1"
#   } else{
#     preds_ensembl[i,6] <- "0"
#   }
#   preds_ensembl[i,7] <- as.character(y[,1][which(y$names==preds_ensembl[i,1])])
#   };rm(i)

# preds_ensembl[,7] <- factor(preds_ensembl[,7],levels = c("1","0"))
# preds_ensembl[,6] <- factor(preds_ensembl[,6],levels = c("1","0"))
# preds_ensembl[,2] <- factor(preds_ensembl[,2],levels = c("1","0"))

preds_df <-   as.data.frame(preds_list2[[1]]$genes)
preds_df$prob_lda <- preds_list2[[1]]$prob_lda
preds_df$prob_svm <- preds_list2[[2]]$prob_svm
preds_df$prob_logit <- preds_list2[[3]]$prob_logit
preds_df$pred_lda <- preds_list2[[1]]$pred_lda
preds_df$pred_svm <- preds_list2[[2]]$pred_svm
preds_df$pred_logit<- preds_list2[[3]]$pred_logit
################################################################################
write.table(preds_df,paste0(dir,"/","ADASAMPLING_WEIGHTED.TSV"),na = "",quote = F,row.names = F,sep = "\t")
################################################################################
# #caret::confusionMatrix(data = factor(preds_ensembl[,6],levels = c("1","0")), reference = preds_ensembl[,7],mode = "everything", positive="1")
# cm_2 <- caret::confusionMatrix(data = factor(preds_ensembl[,2],levels = c("1","0")), reference = preds_ensembl[,7],mode = "everything", positive="1")
# print(cm_2)
# 
# can_sum_pred <- sum(as.numeric(as.character(preds_ensembl[,2])))
# can_sum_obs <- sum(as.numeric(as.character(preds_ensembl[,7])))
# 
# can_sum_pred
# can_sum_obs
# # cm_6 <- caret::confusionMatrix(data = factor(preds_ensembl[,6],levels = c("1","0")), reference = preds_ensembl[,7],mode = "everything", positive="1")
# # print(cm_6)
# 
# ################################################################################
# #https://blog.revolutionanalytics.com/2016/03/com_class_eval_metrics_r.html#perclass
# n = sum(cm_2$table) # number of instances
# nc = nrow(cm_2$table) # number of classes
# diag = diag(cm_2$table) # number of correctly classified instances per class 
# rowsums = apply(cm_2$table, 1, sum) # number of instances per class
# colsums = apply(cm_2$table, 2, sum) # number of predictions per class
# p = rowsums / n # distribution of instances over the actual classes
# q = colsums / n # distribution of instances over the predicted classes
# accuracy = sum(diag) / n 
# precision = diag / colsums 
# recall = diag / rowsums 
# f1 = 2 * precision * recall / (precision + recall) 
# data.frame(precision, recall, f1) 
# 
# macroPrecision = mean(precision)
# macroRecall = mean(recall)
# macroF1 = mean(f1)
# data.frame(macroPrecision, macroRecall, macroF1)
# 
# avgAccuracy = sum(diag(cm_2$table)) / sum(cm_2$table)
# 
# mcIndex = which(rowsums==max(rowsums))[1] # majority-class index
# mcAccuracy = as.numeric(p[mcIndex]) 
# mcRecall = 0*p;  mcRecall[mcIndex] = 1
# mcPrecision = 0*p; mcPrecision[mcIndex] = p[mcIndex]
# mcF1 = 0*p; mcF1[mcIndex] = 2 * mcPrecision[mcIndex] / (mcPrecision[mcIndex] + 1)
# data.frame(mcRecall, mcPrecision, mcF1) 
# 
# expAccuracy = sum(p*q)
# kappa = (accuracy - expAccuracy) / (1 - expAccuracy)
# kappa
# ################################################################################
# colnames(preds_ensembl) <- c("genes","LDA","SVM","LOGIT","NA","ENSEMBLE","CANDIDATES")
#write.table(preds_ensembl,paste0(dir,"/","ADASAMPLING.TSV"),na = "",quote = F,row.names = F,sep = "\t")
# 
# # require(PRROC)
# # 
# # 
# # 
# # roc<-roc.curve(scores.class0 = preds[,1], scores.class1 = preds[,2],curve = T)
# # pr<-pr.curve(scores.class0 = preds[,1], scores.class1 = preds[,2],curve = T)
# # 
# # plot(pr)
# # plot(roc)
# 
# 
# 
# ################################################################################
# # preds_ensembl$final_class <- NA
# # for(i in 1:nrow(preds_ensembl)){
# #   if(is.na(preds_ensembl[1,2])){
# #     preds_ensembl[1,7] <- N
# #   }
# #   
# # }
# 
# factoextra::fviz_pca_ind(PCC,
#              label = "none", # hide individual labels
#              habillage =preds_ensembl[,2] , # color by groups
#              palette = c("#00AFBB", "#FC4E07"),
#              addEllipses = F # Concentration ellipses
# )
# 
# dat <- data.frame(obs=preds_ensembl[,7],
#                   pred=preds_ensembl[,2],
#                   class1= preds_list[[1]][,2],
#                   class0= preds_list[[1]][,1]
# )
# colnames(dat)[3:4] <- c("1","0")
# defaultSummary(dat, lev = c("1","0"))
# twoClassSummary(dat, lev = c("1","0"))
# prSummary(dat,lev = c("1","0"))
# mnLogLoss(dat,lev = c("1","0"))
# 
