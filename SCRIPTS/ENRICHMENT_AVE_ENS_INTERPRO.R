require(readxl)
# require(ggplot2)
require(riceidconverter)
# require(gprofiler2)
# require(openxlsx)
# require(UpSetR)
# require(ggpubr)
require(bc3net)


bc3net_enrichment_interpro <- function(dir,filename,sheet,sep,field,out_dir,out_name,file_interpro){
  
  
file_interpro <- read.table(file_interpro,sep = sep,header = T)

unique_interpro <- unique(file_interpro$motif_id)

unique_interpro_list <- list()

for(i in 1:length(unique_interpro)){
  unique_interpro_list[[i]] <- unique(file_interpro$gene_id[which(file_interpro$motif_id==unique_interpro[[i]])])
};rm(i)

names(unique_interpro_list) <- unique_interpro



    WEIGHTED_DF <-
      as.data.frame(read_xlsx(paste0(dir, "/", filename),sheet = sheet))

    universe = riceidconverter::RiceIDConvert(WEIGHTED_DF$proteins, fromType = "MSU", toType = "RAP")
    universe = unique(universe$RAP)
    universe = universe[which(universe != "None")]
    universe = as.character(universe)
    
    
    WEIGHTED_DF_datasets <- unique(WEIGHTED_DF[,field])
    WEIGHTED_DF_datasets <- WEIGHTED_DF_datasets[which(WEIGHTED_DF_datasets!="N/A")]
    results_list <- list()
    
    
    message("TO PROCESS...")
    #print(WEIGHTED_DF_datasets)
    print(tapply(WEIGHTED_DF[,field], WEIGHTED_DF[,field], length))
    message("...")
    
    tab.hypg_list <- list()
    
    
    for(j in 1:length(WEIGHTED_DF_datasets)){
      message(WEIGHTED_DF_datasets[[j]])
      #j <- 2
      message(j)
      #j <- 1
      candidates_gene_ids <-
        riceidconverter::RiceIDConvert(WEIGHTED_DF$proteins[which(WEIGHTED_DF[,field]==WEIGHTED_DF_datasets[[j]])], fromType = "MSU", toType = "RAP")
      candidates_gene_ids <- unique(candidates_gene_ids$RAP)
      candidates_gene_ids <-
        candidates_gene_ids[which(candidates_gene_ids != "None")]
      candidates_gene_ids <- as.character(candidates_gene_ids)
      
      tab.hypg <- enrichment(candidates_gene_ids, universe, unique_interpro_list, verbose=F,adj = "fdr")
      tab.hypg$approach <- WEIGHTED_DF_datasets[j]
      tab.hypg <- tab.hypg[which(tab.hypg$padj<0.05),]
      
      tab.hypg_list[[j]] <- tab.hypg
      
    };rm(j)
      
    tab.hypg_list <- tab.hypg_list[!sapply(tab.hypg_list,is.null)]
    tab.hypg_list <- do.call(rbind,tab.hypg_list)
    tab.hypg_list$description <- NA
    
    unique_ids <- unique(tab.hypg_list$TermID)
    
    for(i in 1:length(unique_ids)){
      x_f <- file_interpro[which(file_interpro$motif_id==unique_ids[i]),]
      tab.hypg_list$description[tab.hypg_list$TermID %in% unique_ids[i]] <- paste(unique(x_f$description),collapse = "//")
      
    };rm(i)
    
    openxlsx::write.xlsx(tab.hypg_list, file = paste0(out_dir,"/",out_name,".xlsx"))
    message("DONE!")
    return(tab.hypg_list)
    
    }

    #tab.hypg_list$dataset <- field

################################################################################
file_interpro <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS/interpro.osa.tsv"
dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
out_dir <- dir
filename = "COMBINED_RESULTS.xlsx"
sheet = "Candidate ranks"
field = "status D1"
out_name <- "GO_GOCOMPARE_RANKING_D1_INTERPRO"
sep <- "\t"
x1 <- bc3net_enrichment_interpro(dir,filename,sheet,sep,field,out_dir,out_name,file_interpro)
################################################################################
file_interpro <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS/interpro.osa.tsv"
dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
out_dir <- dir
filename = "COMBINED_RESULTS.xlsx"
sheet = "Candidate ranks"
field = "status D2"
out_name <- "GO_GOCOMPARE_RANKING_D2_INTERPRO"
sep <- "\t"
x2 <- bc3net_enrichment_interpro(dir,filename,sheet,sep,field,out_dir,out_name,file_interpro)
################################################################################
file_interpro <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS/All-pathways-of-O.-sativa-Japonica-Group.csv"
dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
out_dir <- dir
filename = "COMBINED_RESULTS.xlsx"
sheet = "Candidate ranks"
field = "status D1"
out_name <- "GO_GOCOMPARE_RANKING_D1_ORYZACYC"
sep <- "|"
x11 <- bc3net_enrichment_interpro(dir,filename,sheet,sep,field,out_dir,out_name,file_interpro)
################################################################################
file_interpro <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/INPUTS/All-pathways-of-O.-sativa-Japonica-Group.csv"
dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
out_dir <- dir
filename = "COMBINED_RESULTS.xlsx"
sheet = "Candidate ranks"
field = "status D2"
out_name <- "GO_GOCOMPARE_RANKING_D2_ORYZACYC"
sep <- "|"
x22 <- bc3net_enrichment_interpro(dir,filename,sheet,sep,field,out_dir,out_name,file_interpro)
################################################################################