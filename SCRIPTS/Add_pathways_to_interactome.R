require(readxl);require(riceidconverter);require(openxlsx)

file_interpro <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS/All-pathways-of-O.-sativa-Japonica-Group.csv"
dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
out_dir <- dir
filename = "COMBINED_RESULTS.xlsx"
sheet = "Candidate ranks"
out_name <- "COMBINED_RESULTS_PATHWAYS"
out_name2 <- "COMBINED_RESULTS_PATHWAYS_IND"
file_interpro <- read.table(file_interpro,sep = "|",header = T)

WEIGHTED_DF <-
  as.data.frame(read_xlsx(paste0(dir, "/", filename),sheet = sheet))

WEIGHTED_DF$pathways <- NA

WEIGHTED_DF_list <- list()
for(i in 1:nrow(WEIGHTED_DF)){
 # i <- 74
  #i <- 11190
  message(i)
  x <- riceidconverter::RiceIDConvert(WEIGHTED_DF$proteins[[i]], fromType = "MSU", toType = "RAP")
  x_status <- file_interpro$description[file_interpro$gene_id %in% as.character(x$RAP)]
  if(length(x_status)>0){
    message(paste(i," ..."))
    x_i  <- WEIGHTED_DF[rep(i, length(x_status)),]
    x_i$pathways <- unique(x_status)
    WEIGHTED_DF_list[[i]] <- x_i
    #WEIGHTED_DF_list[[i]] <- 
    WEIGHTED_DF$pathways[[i]] <- paste(unique(x_status),collapse = "//")    
  } else {
    WEIGHTED_DF$pathways[[i]] <- NA 
    WEIGHTED_DF_list[[i]] <- WEIGHTED_DF[i,]
  }
  #x_status <- WEIGHTED_DF$`status D1`[WEIGHTED_DF$proteins %in% x$MSU]
};rm(i)

WEIGHTED_DF_list <- do.call(rbind,WEIGHTED_DF_list)


openxlsx::write.xlsx(WEIGHTED_DF, file = paste0(out_dir,"/",out_name,".xlsx"))
openxlsx::write.xlsx(WEIGHTED_DF_list, file = paste0(out_dir,"/",out_name2,".xlsx"))
# write.table(WEIGHTED_DF, file = paste0(out_dir,"/",out_name,".tsv"),sep = "\t",row.names = F)
# write.table(WEIGHTED_DF_list, file = paste0(out_dir,"/",out_name2,".tsv"),sep = "\t",row.names = F)
