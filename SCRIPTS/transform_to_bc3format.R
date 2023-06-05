require(readxl)
file_interpro <- as.data.frame(readxl::read_xlsx("D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS/All-pathways-of-O.-sativa-Japonica-Group.xlsx"))

pathways <- list()
for(i in 1:nrow(file_interpro)){
  #i <- 1
  message(i)
  x <- file_interpro$`Genes of pathway`[[i]]
  x <- strsplit(x,split = " // ")[[1]]
  x <-x[x!=""]
  
  if(length(x)>0){
  pathways[[i]] <- data.frame(gene_id=x,motif_id=file_interpro$Pathway_id[[i]],description=file_interpro$Pathways[[i]])
  } else {
    pathways[[i]] <- data.frame(gene_id=NA,motif_id=file_interpro$Pathway_id[[i]],description=file_interpro$Pathways[[i]])
    
  }
};rm(i)

pathways <- do.call(rbind,pathways)
pathways <- pathways[which(!is.na(pathways$gene_id)),]
pathways <- pathways[which(pathways$gene_id!=""),]
write.table(pathways,"D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS/All-pathways-of-O.-sativa-Japonica-Group.csv",
            na = "",
            quote = T,sep = "|",row.names = F
            )
