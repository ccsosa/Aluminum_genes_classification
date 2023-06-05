require(readxl);require(ggplot2);require(agricolae);require(reshape2);require(ggpubr)
dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
################################################################################
data_LFC1 <- as.data.frame(readxl::read_xlsx(paste0(dir,"/GOCOMPARE/","CANDIDATES_GOCompare_ENSEMBLE_UNWEIGHTED_ALL.xlsx"),sheet = "UNWEIGHTED"))
data_LFC1 <- data_LFC1[which(data_LFC1$status==0),]
row.names(data_LFC1) <- data_LFC1$genes;data_LFC1$genes <- NULL

data_LFC1_W <- as.data.frame(readxl::read_xlsx(paste0(dir,"/GOCOMPARE/","CANDIDATES_GOCompare_ENSEMBLE_WEIGHTED_ALL.xlsx"),sheet = "WEIGHTED"))
data_LFC1_W <- data_LFC1[which(data_LFC1_W$status==0),]
row.names(data_LFC1_W) <- data_LFC1_W$genes;data_LFC1_W$genes <- NULL


dd <- colnames(data_LFC1)
prob_cols <- dd[c(9:11)]
pred_cols <- dd[c(12:14)]

data_LFC1_list <- list()
data_LFC1_W_list <- list()

for(i in 1:length(prob_cols)){
    data_LFC1_list[[i]] <- data.frame(
    variable= prob_cols[i],
      value= data_LFC1[,prob_cols[i]][which(data_LFC1[,pred_cols[i]]==1)],
      approach= "Unweighted")
    
    data_LFC1_W_list[[i]] <- data.frame(
      variable= prob_cols[i],
      value= data_LFC1_W[,prob_cols[i]][which(data_LFC1_W[,pred_cols[i]]==1)],
      approach= "Weighted")
    
};rm(i)
data_LFC1 <- do.call(rbind,data_LFC1_list)
data_LFC1_W <- do.call(rbind,data_LFC1_W_list)


data <- rbind(data_LFC1,data_LFC1_W)

data$ap_meth <- paste0(data$variable,"-",data$approach)
data$approach <- factor(data$approach,levels = c("Unweighted","Weighted"))
# comparison<-agricolae::kruskal(data$value,data$ap_meth,group=F, main="data_LFC1",alpha = 0.05,p.adj = "BH")

################################################################################
data_LFC2 <- as.data.frame(readxl::read_xlsx(paste0(dir,"/JJG_DEG_1/","CANDIDATES_LFC1_ENSEMBLE_UNWEIGHTED_ALL.xlsx"),sheet = "UNWEIGHTED"))
data_LFC2 <- data_LFC2[which(data_LFC2$status==0),]
row.names(data_LFC2) <- data_LFC2$genes;data_LFC2$genes <- NULL

data_LFC2_W <- as.data.frame(readxl::read_xlsx(paste0(dir,"/JJG_DEG_1/","CANDIDATES_LFC1_ENSEMBLE_WEIGHTED_ALL.xlsx"),sheet = "WEIGHTED"))
data_LFC2_W <- data_LFC2_W[which(data_LFC2_W$status==0),]
row.names(data_LFC2_W) <- data_LFC2_W$genes;data_LFC2_W$genes <- NULL

dd2 <- colnames(data_LFC2)
prob_cols <- dd2[c(9:11)]
pred_cols <- dd2[c(12:14)]

data_LFC2_list <- list()
data_LFC2_W_list <- list()

for(i in 1:length(prob_cols)){
  data_LFC2_list[[i]] <- data.frame(
    variable= prob_cols[i],
    value= data_LFC2[,prob_cols[i]][which(data_LFC2[,pred_cols[i]]==1)],
    approach= "Unweighted")
  
  data_LFC2_W_list[[i]] <- data.frame(
    variable= prob_cols[i],
    value= data_LFC2_W[,prob_cols[i]][which(data_LFC2_W[,pred_cols[i]]==1)],
    approach= "Weighted")
  
};rm(i)
data_LFC2 <- do.call(rbind,data_LFC2_list)
data_LFC2_W <- do.call(rbind,data_LFC2_W_list)


#data_LFC1_W$variable <- paste(data_LFC1_W$variable," (WEIGHTED)")

data2 <- rbind(data_LFC2,data_LFC2_W)

data2$ap_meth <- paste0(data2$variable,"-",data2$approach)
data2$approach <- factor(data2$approach,levels = c("Unweighted","Weighted"))

#unique(data2$ap_meth)
################################################################################
#un_ad <- unique(data$ap_meth)

#comparison <- 

data <- data[ data$variable %in% dd[c(9:11)],]

# #p <- ggboxplot

data2 <- data2[ data2$variable %in% dd2[c(9:11)],]

# data$ap_meth <- factor(data$ap_meth,levels = c("PROM_PROB_PULEARN-Unweighted",
#                                                "PROM_PROB_PULEARN-Weighted",
#                                                "PROM_PROB_ADASAMPLING-Unweighted",
#                                                "PROM_PROB_ADASAMPLING-Weighted",
#                                                "PROM_PROB_ALL-Unweighted",
#                                                "PROM_PROB_ALL-Weighted" ))

# data2$ap_meth <- factor(data2$ap_meth,levels = c("PROM_PROB_PULEARN-Unweighted",
#                                                "PROM_PROB_PULEARN-Weighted",
#                                                "PROM_PROB_ADASAMPLING-Unweighted",
#                                                "PROM_PROB_ADASAMPLING-Weighted",
#                                                "PROM_PROB_ALL-Unweighted",
#                                                "PROM_PROB_ALL-Weighted" ))

################################################################################

dict <- data.frame(from=c("PROM_PROB_PULEARN-Unweighted",
                          "PROM_PROB_PULEARN-Weighted",
                          "PROM_PROB_ADASAMPLING-Unweighted",
                          "PROM_PROB_ADASAMPLING-Weighted",
                          "PROM_PROB_ALL-Unweighted",
                          "PROM_PROB_ALL-Weighted" ),
                   to =c("PULearning-Unweighted",
                         "PULearning-Weighted",
                         "AdaSampling-Unweighted",
                         "AdaSampling-Weighted",
                         "Average-Unweighted",
                         "Average-Weighted" ))




for(i in 1:nrow(dict)){
  data$ap_meth[data$ap_meth %in% dict$from[i]] <- dict$to[[i]]
  data2$ap_meth[data2$ap_meth %in% dict$from[i]] <- dict$to[[i]]
};rm(i)

#x <- combn(unique(dd[c(9:11)]),m = 2,simplify = T)
x <- combn(c("PULearning-Unweighted",
             "PULearning-Weighted",
             "AdaSampling-Unweighted",
             "AdaSampling-Weighted",
             "Average-Unweighted",
             "Average-Weighted" ),m = 2,simplify = T)

x <- x[,c(1,10,15)]
list_x <- list()
for(i in 1:ncol(x)){
  list_x[[i]] <- c(as.character(x[1,i]),as.character(x[2,i]))
}
################################################################################
data$ap_meth <- factor(data$ap_meth,levels = c("PULearning-Unweighted",
                                               "PULearning-Weighted",
                                               "AdaSampling-Unweighted",
                                               "AdaSampling-Weighted",
                                               "Average-Unweighted",
                                               "Average-Weighted" ))

data2$ap_meth <- factor(data2$ap_meth,levels = c("PULearning-Unweighted",
                                                 "PULearning-Weighted",
                                                 "AdaSampling-Unweighted",
                                                 "AdaSampling-Weighted",
                                                 "Average-Unweighted",
                                                 "Average-Weighted" ))
################################################################################
p1 <- ggviolin(data, x = "ap_meth", y = "value",fill="approach",
               #color = "approach",
               add = c("boxplot","mean"),
               palette = c("green","purple"),
               add.params = list(fill = "white"),
               xlab = "",
               ylab = "Probability",legend.title="",
               font.caption = c(20, "bold"))

p1 <- p1 + stat_compare_means(comparisons = list_x,size=9)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.1,size=9)  +
  rotate_x_text(90)# Add global p-value
#ggpar(p, ylim = c(0,1))
#p1 <- ggpar(p1,ylim = c(0,1))
p1 <- p1 + rremove("x.text")+
  theme(text = element_text(size = 30))

p2<- ggviolin(data2, x = "ap_meth", y = "value",fill="approach",
              #color = "approach",
              add = c("boxplot","mean"),
              palette = c("green","purple"),
              add.params = list(fill = "white"),
              xlab = "",
              ylab = "Probability",legend.title="",
              font.caption = c(20, "bold"))

p2 <- p2 + stat_compare_means(comparisons = list_x,size=9)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 0.1,size=9)  +
  rotate_x_text(90)# Add global p-value
#ggpar(p, ylim = c(0,1))
#p2 <- ggpar(p2,ylim = c(0,1))+
p2 <- p2 +
  theme(text = element_text(size = 30),
  ) # change text size of theme components


p_all <-ggarrange(
  p1, p2, labels = c("A", "B"),ncol =1,nrow = 2,
  common.legend = TRUE, legend = "bottom",font.label = c(size=40)
)

ggexport(p_all,
         filename = paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS","/","PROB_VIOLIN_CAND.pdf"),
         width = 18,
         height = 25,
         #pointsize = 400,
         res =1000,
         verbose = TRUE
)
#######################################