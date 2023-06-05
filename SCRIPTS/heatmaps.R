require(readxl);require(ggplot2);library(vegan);library(ggheatmap);library(pheatmap)
################################################################################
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
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
dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"

data <- as.data.frame(readxl::read_xlsx(paste0(dir,"/","COMBINED_RESULTS.xlsx")))
row.names(data) <- data$proteins;genes <- data$proteins
data$proteins <- NULL

#colnames(data)
data <- data[,c(1,3,4,2,5,6)]
# data <- t(data)
# colnames(data) <- genes
# xx <- vegdist(data, method = "bray")
# plot(hclust(xx,method = "ward.D"))

mat_num <- matrix(ncol = 6,nrow = 6)
colnames(mat_num) <- colnames(data)
row.names(mat_num) <- colnames(data)
for(i in 1:nrow(mat_num)){
  for(j in 1:ncol(mat_num)){
    mat_num[i,j] <- jaccard(row.names(data[which(data[,i]==1),]),row.names(data[which(data[,j]==1),]))
  };rm(j)
};rm(i)
#heatmap(as.matrix(mat_num))

row.names(mat_num)

# row_metaData <- data.frame(dataset=c(base::rep("dataset 1",3),
#                                      base::rep("dataset 2",3)))

 row_metaData <- data.frame(dataset=c("dataset 1","dataset 1","dataset 1",
                                      "dataset 2","dataset 2","dataset 2"))
rownames(row_metaData) <- colnames(mat_num)



# exprcol <- c("#EE7E30","#5D9AD3")
# names(exprcol) <- c("dataset 1","dataset 2")
# col <- list(dataset=exprcol)

mat_num[lower.tri(mat_num)] <- NA

mat_num <- as.data.frame(mat_num)
mat_num2 <- round(mat_num,3)
mat_num2[is.na(mat_num2)] <- ""
#mat_num2 <- as.character(mat_num2)
diag(mat_num) <- NA
diag(mat_num2) <- ""


p <- pheatmap(mat_num,cluster_rows = F,cluster_cols = F,clustering_distance_rows = "euclidean",clustering_method = "ward.D",
              annotation_colors = list(dataset=c("dataset 1"="red","dataset 2"="blue")),
              annotation_col=row_metaData,legend_labels = T,display_numbers = mat_num2,border_color = "gray",fontsize_number = 12,#)#,
               color=colorRampPalette(c("purple","violet", "white", "yellowgreen"))(5))
              #annotation_row =  row_metaData,
              #annotation_color = col,
              #border = "grey"#,legendName = "Jaccard index"
#              annotation_color = exprcol
p
save_pheatmap_pdf(x = p,
                  filename =paste0("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS","/","HEATMAP_JACCARD.pdf"),
                  width = 8,
                  height = 8 )
