require(UpSetR);require(xlsx);require(readxl)


data <- as.data.frame(read_excel("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE/COMBINED_RESULTS.xlsx"))
row.names(data) <- data$proteins
data$proteins <- NULL
upset(data,sets = colnames(data), mb.ratio = c(0.55, 0.45), order.by = "freq", 
      keep.order = TRUE,
      #group.by = "sets",
      sets.x.label = "Genes",set_size.show = TRUE,decreasing = T,set_size.numbers_size = T)
      


x_up <- upset(data,sets = colnames(data), mb.ratio = c(0.55, 0.45), order.by = "freq", 
      keep.order = TRUE,
      #group.by = "sets",
      sets.x.label = "Proteins",set_size.show = TRUE,decreasing = T,set_size.numbers_size = T,
      nintersects = 20,line.size = 1.3,point.size = 4,text.scale = 1.5)


dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS"
#saving upset plot in pdf format
pdf(paste0(dir,"/","UPSET_UP.pdf"),height = 8, width = 12)#15
x_up
dev.off()


##dataset1

data1 <- data[,c("dataset 1 (728 proteins)","Unweighted graph predicted proteins (dataset 1)","Weighted graph predicted proteins (dataset 1)")]
x_up <- upset(data1,sets = colnames(data1), mb.ratio = c(0.55, 0.45), order.by = "freq", 
              keep.order = TRUE,
              #group.by = "sets",
              sets.x.label = "Proteins",set_size.show = TRUE,decreasing = T,set_size.numbers_size = T,
              line.size = 1.3,point.size = 4,text.scale = 1.5)

dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS"
#saving upset plot in pdf format
pdf(paste0(dir,"/","UPSET_D1.pdf"),height = 8, width = 12)#15
x_up
dev.off()

##dataset2

data2 <- data[,c("dataset 2 (133 proteins)","Unweighted graph predicted proteins (dataset 2)","Weighted graph predicted proteins (dataset 2)")]
x_up <- upset(data2,sets = colnames(data2), mb.ratio = c(0.55, 0.45), order.by = "freq", 
              keep.order = TRUE,
              #group.by = "sets",
              sets.x.label = "Proteins",set_size.show = TRUE,decreasing = T,set_size.numbers_size = T,
              line.size = 1.3,point.size = 4,text.scale = 1.5)

dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS"
#saving upset plot in pdf format
pdf(paste0(dir,"/","UPSET_D2.pdf"),height = 8, width = 12)#15
x_up
dev.off()
# require(riceidconverter)
# 
# genes <- read.table("clipboard")
# ri_df <- riceidconverter::RiceIDConvert(myID = genes,fromType = "MSU",toType = "RAP")
# 
# write.table(ri_df,"D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS/MSU_TO_RAP.tsv",sep = "\t",row.names = F)
