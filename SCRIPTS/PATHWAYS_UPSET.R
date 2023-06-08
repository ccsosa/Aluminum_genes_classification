require(readxl)
require(ggplot2)
require(riceidconverter)
require(gprofiler2)
require(openxlsx)
require(UpSetR)
library(pheatmap)
################################################################################
dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
################################################################################
#function to save pheatmap
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
################################################################################
################################################################################
combined_results <- as.data.frame(
  readxl::read_xlsx(paste0(dir,"/","COMBINED_RESULTS_PATHWAYS_IND.xlsx")))
#row.names(combined_results) <- combined_results$proteins;combined_results$proteins; <- NULL
################################################################################
DATASETS <- c("HIGH CONFIDENCE",
              "TRAINING",
              "MEDIUM CONFIDENCE (Weighted)",
              "MEDIUM CONFIDENCE (Unweighted)",
              "TRAINING (UNWEIGHTED)",
              "TRAINING (WEIGHTED)"
)


dict <- data.frame(from=DATASETS,
                   to= c("High confidence",
                         "Positive class",
                         "Medium confidence (Weighted)",
                         "Medium confidence (Unweighted)",
                         "Positive class (Unweighted)",
                         "Positive class (Weighted)"
                   ))

for(i in 1:nrow(dict)){
  combined_results$`status D1`[combined_results$`status D1` %in% dict$from[i]] <- dict$to[[i]]
  combined_results$`status D2`[combined_results$`status D2` %in% dict$from[i]] <- dict$to[[i]]
  
};rm(i)
################################################################################
combined_results1 <- combined_results[,c("proteins","status D1","pathways")]
combined_results1$dataset <- "Dataset 1"
combined_results2 <- combined_results[,c("proteins","status D2","pathways")]
combined_results2$dataset <- "Dataset 2"

colnames(combined_results1)[2] <- "status";colnames(combined_results2)[2] <- "status"
comb_results_final <- rbind(combined_results1,combined_results2)
comb_results_final <- comb_results_final[which(comb_results_final$status!="N/A"),]
#comb_results_final <- comb_results_final[which(!is.na(comb_results_final$pathways)),]
comb_results_final$label <- paste(comb_results_final$status,"-",comb_results_final$dataset)
################################################################################
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
################################################################################
labels <- unique(comb_results_final$label)
pathways <- unique(comb_results_final$pathways)

mat_num <- as.data.frame(matrix(nrow = length(pathways),ncol= length(labels)))
colnames(mat_num) <- labels
row.names(mat_num) <- c("No available",pathways[-1])


for(i in 1:ncol(mat_num)){
  for(j in 1:nrow(mat_num)){
    if(j ==1){
      mat_num[j,i] <- nrow(comb_results_final[which(is.na(comb_results_final$pathways) &
                               comb_results_final$label==labels[[i]]),])

    } else {
      mat_num[j,i] <- nrow(comb_results_final[which(comb_results_final$pathways==pathways[[j]] &
                               comb_results_final$label==labels[[i]]),])   

    }
  };rm(j)     
  
};rm(i)

#Qcolnames(mat_num)[1] <- "No available"
mat_num <- mat_num[-1,]
mat_num <- mat_num[,-11]
#mat_num <- scale(mat_num,center = T,scale = T)
row_metaData <- data.frame(dataset=c("dataset 1","dataset 1",
                                     "dataset 1","dataset 1",
                                     "dataset 1","dataset 1",
                                     "dataset 2","dataset 2",
                                     "dataset 2","dataset 2",
                                     "dataset 2"))
rownames(row_metaData) <- colnames(mat_num)

mat_num <- mat_num[which(rowSums(mat_num)>0),]
#mat_num[mat_num==0] <- NA
#mat_num[is.na(mat_num)] <- as.double("NA")

#mat_num <- as.data.frame(t(mat_num))
cols <- rev(hcl.colors(24, "Greens"))
cols<- c("gray",cols)

################################################################################
p  <- pheatmap(mat_num,cluster_rows = F,
               cluster_cols = T,
               #clustering_distance_rows = "correlation",
               #clustering_method = "ward.D",
               annotation_colors = list(dataset=c("dataset 1"="red","dataset 2"="blue")),
               annotation_col=row_metaData,
               na_col="gray",
            #legend_labels = T,
            border_color = "black",
            #scale ="column",
            fontsize_number = 18,fontsize_row = 10,fontsize_col = 12,
            #color=colorRampPalette(c("purple","violet", "white", "yellowgreen"))(5),
            color =cols#,
           #legend_breaks = c(0,1,2,4,8,12,16,20), 
           #legend_labels = c("0","1","2","4","8","12","16","20")
           )#,
            #scale = "column",)

save_pheatmap_pdf(x = p,
                  filename =paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS","/","HEATMAP_PATHWAYS.pdf"),
                  width = 12,
                  height = 60)

################################################################################
# mat_num2 <- mat_num#[mat_num==0] <- NA
# mat_num2[mat_num2==0] <- NA

x_p95 <- list()

x_p95 <- lapply(1:ncol(mat_num2),function(i){
  message(i)
  x_i <- mat_num2[which(mat_num2[,i] > quantile(mat_num2[,i],probs = 0.95,na.rm = T)),]
  if(nrow(x_i)>0){
  x_i <- x_i[order(x_i[,i],decreasing = T),]
  x_i <- data.frame(pathway=row.names(x_i),count=x_i[,i],dataset=colnames(x_i)[i])
  return(x_i)
  }
  })

x_p95 <- do.call(rbind,x_p95)
write.csv(x_p95, paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS","/","top95.csv"),na = "",row.names = F)
