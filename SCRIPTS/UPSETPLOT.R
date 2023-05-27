require(UpSetR);require(xlsx);require(readxl)


data <- as.data.frame(read_excel("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE/COMBINED_RESULTS.xlsx"))
row.names(data) <- data$genes
data$genes <- NULL
upset(data,sets = colnames(data), mb.ratio = c(0.55, 0.45), order.by = "freq", 
      keep.order = TRUE,
      #group.by = "sets",
      sets.x.label = "Genes",set_size.show = TRUE,decreasing = T,set_size.numbers_size = T)
      

# require(riceidconverter)
# 
# genes <- read.table("clipboard")
# ri_df <- riceidconverter::RiceIDConvert(myID = genes,fromType = "MSU",toType = "RAP")
# 
# write.table(ri_df,"D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS/MSU_TO_RAP.tsv",sep = "\t",row.names = F)
