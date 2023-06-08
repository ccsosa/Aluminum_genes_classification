require(readxl);require(agricolae);require(ggpubr)

x <- as.data.frame(read.csv("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE/NETWORK_NODE.csv",header = T,row.names = 1))


x_k1 <- kruskal(as.numeric(x$BetweennessCentrality),as.factor(x$status.D1),group = F,p.adj = "fdr")

x_k1

################################################################################
DATASETS <- c("HIGH CONFIDENCE",
              "TRAINING",
              "MEDIUM CONFIDENCE (Weighted)",
              "MEDIUM CONFIDENCE (Unweighted)",
              "TRAINING (UNWEIGHTED)",
              "TRAINING (WEIGHTED)",
              "N/A"
)
################################################################################

dict <- data.frame(from=DATASETS,
                   to= c("High confidence",
                         "Positive class",
                         "Medium confidence (Weighted)",
                         "Medium confidence (Unweighted)",
                         "Positive class (Unweighted)",
                         "Positive class (Weighted)",
                         "Unlabeled"
                   ))

for(i in 1:nrow(dict)){
  x$status.D1[x$status.D1 %in% dict$from[i]] <- dict$to[[i]]
  x$status.D2[x$status.D2 %in% dict$from[i]] <- dict$to[[i]]
  
};rm(i)
################################################################################
x$status.D1 <- factor(x$status.D1,levels = c("Positive class",
                                             "Positive class (Unweighted)",
                                             "Positive class (Weighted)",
                                             "High confidence",
                                             "Medium confidence (Weighted)",
                                             "Medium confidence (Unweighted)",
                                             "Unlabeled"))

x$status.D2 <- factor(x$status.D2,levels = c("Positive class",
                                             "Positive class (Unweighted)",
                                             "Positive class (Weighted)",
                                             "High confidence",
                                             "Medium confidence (Weighted)",
                                             "Medium confidence (Unweighted)",
                                             "Unlabeled"))
################################################################################
#x <- combn(unique(dd[c(9:11)]),m = 2,simplify = T)
x_comb <- combn(unique(x$status.D1),m = 2,simplify = T)
x_comb <- x_comb[,c(1:6)]
#x <- x[,c(1,10,15)]
list_x <- list()
for(i in 1:ncol(x_comb)){
  list_x[[i]] <- c(as.character(x_comb[1,i]),as.character(x_comb[2,i]))
}

################################################################################
################################################################################
p1 <- ggboxplot(x, x = "status.D1", y = "Degree",fill="status.D1",
               #color = "status.D1",
               #add = c("boxplot","mean"),
              # palette = c("green","purple"),
               add.params = list(fill = "white"),
              outlier.shape = NA,
              #trim = T,
               xlab = "",
              order=c("Positive class",
                      "Positive class (Unweighted)",
                      "Positive class (Weighted)",
                      "High confidence",
                      "Medium confidence (Weighted)",
                      "Medium confidence (Unweighted)",
                    "Unlabeled"),
               ylab = "log10 (Degree)",legend.title="",
               font.caption = c(10, "bold"))
              #bxp.errorbar.width = 0.15)#,size = 0.1,font.label = list(size=10))
#p1
p1 <- p1 + stat_compare_means(comparisons = list_x,size=9)+ # Add pairwise comparisons p-value
  stat_compare_means(size=9)  +
  rotate_x_text(90)# Add global p-value
p1 <- 
  ggpar(p1,
        font.x = c(30, "bold", "black"),
        font.y = c(30, "bold", "black"),
        font.xtickslab = c(22, "bold", "black"),
        font.ytickslab = c(22, "bold", "black"),legend = "none",yscale = "log10")
#p1
ggexport(p1,
         filename = paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE","/","VIOLIN_DEGREE_D1.pdf"),
         width = 18,
         height = 25,
         #pointsize = 400,
         res =1000,
         verbose = TRUE
)
################################################################################
p1 <- ggboxplot(x, x = "status.D1", y = "BetweennessCentrality",fill="status.D1",
               #color = "status.D1",
               #add = c("boxplot","mean"),
               # palette = c("green","purple"),
               add.params = list(fill = "white"),
               #trim = T,
               xlab = "",
               outlier.shape = NA,
               order=c("Positive class",
                       "Positive class (Unweighted)",
                       "Positive class (Weighted)",
                       "High confidence",
                       "Medium confidence (Weighted)",
                       "Medium confidence (Unweighted)",
                       "Unlabeled"),
               ylab = "log 10 (Betweenness Centrality)",legend.title="",
               font.caption = c(10, "bold"))#, orientation = "horiz",size = 0.1,font.label = list(size=10))
#p1
p1 <- p1 + stat_compare_means(comparisons = list_x,size=9)+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 20000,size=9)  +
  rotate_x_text(90)# Add global p-value
p1 <- 
  ggpar(p1,
        font.x = c(30, "bold", "black"),
        font.y = c(30, "bold", "black"),
        font.xtickslab = c(22, "bold", "black"),
        font.ytickslab = c(22, "bold", "black"),legend = "none",yscale = "log10")
ggexport(p1,
         filename = paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE","/","VIOLIN_BetweennessCentrality_D1.pdf"),
         width = 18,
         height = 25,
         #pointsize = 400,
         res =1000,
         verbose = TRUE
)

###############################################################################
p1 <- ggboxplot(x, x = "status.D1", y = "ClusteringCoefficient",fill="status.D1",
               #color = "status.D1",
               #add = c("boxplot","mean"),
               # palette = c("green","purple"),
               add.params = list(fill = "white"),
               outlier.shape = NA,
               #trim = T,
               xlab = "",
               order=c("Positive class",
                       "Positive class (Unweighted)",
                       "Positive class (Weighted)",
                       "High confidence",
                       "Medium confidence (Weighted)",
                       "Medium confidence (Unweighted)",
                       "Unlabeled"),
               ylab = "log10 (Clustering Coefficient)",legend.title="",
               font.caption = c(10, "bold"))#, orientation = "horiz",size = 0.1,font.label = list(size=10))
#p1
p1 <- p1 + stat_compare_means(comparisons = list_x,size=9)+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 20000,size=9)  +
  rotate_x_text(90)# Add global p-value
p1 <- 
  ggpar(p1,
        font.x = c(30, "bold", "black"),
        font.y = c(30, "bold", "black"),
        font.xtickslab = c(22, "bold", "black"),
        font.ytickslab = c(22, "bold", "black"),legend = "none",yscale = "log10")
ggexport(p1,
         filename = paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE","/","VIOLIN_ClusteringCoefficient_D1.pdf"),
         width = 18,
         height = 25,
         #pointsize = 400,
         res =1000,
         verbose = TRUE
)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
p1 <- ggboxplot(x, x = "status.D2", y = "Degree",fill="status.D2",
               #color = "status.D2",
               #add = c("boxplot","mean"),
               # palette = c("green","purple"),
               add.params = list(fill = "white"),
               #trim = T,
               xlab = "",
               order=c("Positive class",
                       "Positive class (Unweighted)",
                       "Positive class (Weighted)",
                       "High confidence",
                       "Medium confidence (Weighted)",
                       "Medium confidence (Unweighted)",
                       "Unlabeled"),
               ylab = "log10 (Degree)",legend.title="",
               font.caption = c(10, "bold"))#, orientation = "horiz",size = 0.1,font.label = list(size=10))
#p1
p1 <- p1 + stat_compare_means(comparisons = list_x,size=9)+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 20000,size=9)  +
  rotate_x_text(90)# Add global p-value
p1 <- 
  ggpar(p1,
        font.x = c(30, "bold", "black"),
        font.y = c(30, "bold", "black"),
        font.xtickslab = c(22, "bold", "black"),
        font.ytickslab = c(22, "bold", "black"),legend = "none",yscale = "log10")
ggexport(p1,
         filename = paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE","/","VIOLIN_DEGREE_D2.pdf"),
         width = 18,
         height = 25,
         #pointsize = 400,
         res =1000,
         verbose = TRUE
)
################################################################################
p1 <- ggboxplot(x, x = "status.D2", y = "BetweennessCentrality",fill="status.D2",
               #color = "status.D2",
               #add = c("boxplot","mean"),
               # palette = c("green","purple"),
               add.params = list(fill = "white"),
              # trim = T,
               xlab = "",
               order=c("Positive class",
                       "Positive class (Unweighted)",
                       "Positive class (Weighted)",
                       "High confidence",
                       "Medium confidence (Weighted)",
                       "Medium confidence (Unweighted)",
                       "Unlabeled"),
               ylab = "log10 (Betweenness Centrality)",legend.title="",
               font.caption = c(10, "bold"))#, orientation = "horiz",size = 0.1,font.label = list(size=10))
#p1
p1 <- p1 + stat_compare_means(comparisons = list_x,size=9)+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 20000,size=9)  +
  rotate_x_text(90)# Add global p-value
p1 <- 
  ggpar(p1,
        font.x = c(30, "bold", "black"),
        font.y = c(30, "bold", "black"),
        font.xtickslab = c(22, "bold", "black"),
        font.ytickslab = c(22, "bold", "black"),legend = "none",yscale = "log10")
ggexport(p1,
         filename = paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE","/","VIOLIN_BetweennessCentrality_D2.pdf"),
         width = 18,
         height = 25,
         #pointsize = 400,
         res =1000,
         verbose = TRUE
)

###############################################################################
p1 <- ggboxplot(x, x = "status.D2", y = "ClusteringCoefficient",fill="status.D2",
               #color = "status.D2",
               #add = c("boxplot","mean"),
               # palette = c("green","purple"),
               add.params = list(fill = "white"),
              # trim = T,
               xlab = "",
               order=c("Positive class",
                       "Positive class (Unweighted)",
                       "Positive class (Weighted)",
                       "High confidence",
                       "Medium confidence (Weighted)",
                       "Medium confidence (Unweighted)",
                       "Unlabeled"),
               ylab = "log10 (Clustering Coefficient)",legend.title="",
               font.caption = c(10, "bold"))#, orientation = "horiz",size = 0.1,font.label = list(size=10))
#p1
p1 <- p1 + stat_compare_means(comparisons = list_x,size=9)+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 20000,size=9)  +
  rotate_x_text(90)# Add global p-value
p1 <- 
  ggpar(p1,
        font.x = c(30, "bold", "black"),
        font.y = c(30, "bold", "black"),
        font.xtickslab = c(22, "bold", "black"),
        font.ytickslab = c(22, "bold", "black"),legend = "none",yscale = "log10")
ggexport(p1,
         filename = paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE","/","VIOLIN_ClusteringCoefficient_D2.pdf"),
         width = 18,
         height = 25,
         #pointsize = 400,
         res =1000,
         verbose = TRUE
)
