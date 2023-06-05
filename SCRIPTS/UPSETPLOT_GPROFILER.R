require(readxl)
require(ggplot2)
require(riceidconverter)
require(gprofiler2)
require(openxlsx)
require(UpSetR)
require(ggpubr)
dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"

file1 <- "GO_GOCOMPARE_RANKING_D1.xlsx"
file2 <- "GO_GOCOMPARE_RANKING_D2.xlsx"

file1 <- as.data.frame(readxl::read_xlsx(paste0(dir,"/",file1)))
file2 <- as.data.frame(readxl::read_xlsx(paste0(dir,"/",file2)))


file1$dataset <- "Dataset 1"
file2$dataset <- "Dataset 2"

funcs <- c("GO:BP","GO:MF","KEGG")

DATASETS <- c("HIGH CONFIDENCE",
              "TRAINING",
              "MEDIUM CONFIDENCE (Weighted)",
              "MEDIUM CONFIDENCE (Unweighted)",
              "TRAINING (UNWEIGHTED)",
              "TRAINING (WEIGHTED)"
)


x_up <- list()
data_func_list <- list()
for(a in 1:length(funcs)){
f_d1 <- file1[which(file1$source==funcs[[a]]),] 
f_d2 <- file2[which(file2$source==funcs[[a]]),] 
f_data <- rbind(f_d1,f_d2)

dict <- data.frame(from=DATASETS,
                   to= c("High confidence",
                         "Positive class",
                         "Medium confidence (Weighted)",
                         "Medium confidence (Unweighted)",
                         "Positive class (Unweighted)",
                         "Positive class (Weighted)"
                   ))

for(i in 1:nrow(dict)){
  f_data$approach[f_data$approach %in% dict$from[i]] <- dict$to[[i]]
};rm(i)

f_data$label <- paste(f_data$approach,"-",f_data$dataset)
term_name <- unique(f_data$term_name)

data_func_mat <- as.data.frame(matrix(ncol = length(unique(f_data$label)),nrow = length(term_name)))
row.names(data_func_mat) <- term_name
colnames(data_func_mat) <- unique(f_data$label)

for(i in 1:nrow(data_func_mat)){
  
  for(j in 1:ncol(data_func_mat)){
    
    x <- f_data[which(f_data$label==colnames(data_func_mat)[j] & f_data$term_name==row.names(data_func_mat)[i]),]
    if(nrow(x)>0){
      data_func_mat[i,j] <- 1
    } else {
      data_func_mat[i,j] <- 0
    }
  };rm(j)
};rm(i)
data_func_list[[a]] <- data_func_mat

x_up[[a]] <- upset(data_func_mat, mb.ratio = c(0.55, 0.45), order.by = "freq", 
              keep.order = TRUE,
              #group.by = "sets",
              sets.x.label = paste("Shared GO term (",funcs[a],")"),set_size.show = TRUE,decreasing = T,set_size.numbers_size = T,
              line.size = 1.3,point.size = 4,text.scale = 1.5)
rm(data_func_mat)

};rm(a)

#ggarrange(x_up[[1]],x_up[[2]],x_up[[3]],ncol = 1,nrow = 3,labels = c("A","B","C"))
x_up[[1]]
x_up[[2]]
x_up[[3]]
