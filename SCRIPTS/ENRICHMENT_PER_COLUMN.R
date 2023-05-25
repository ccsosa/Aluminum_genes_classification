require(readxl)
require(ggplot2)
require(riceidconverter)
require(gprofiler2)
require(openxlsx)


ENRICHMENT_CANDIDATES_FUNC_COL <- function(dir,filename,out_dir,out_name){
  source("https://raw.githubusercontent.com/ccsosa/reducereduntGO/master/reducereduntGO_ALL.R")
  

  WEIGHTED_DF <-
    read_xlsx(paste0(dir, "/", filename))
  
  # training <- WEIGHTED_DF[which(WEIGHTED_DF$status == 1), ]
  # training_gene_ids <-
  #   riceidconverter::RiceIDConvert(training$genes, fromType = "MSU", toType = "RAP")
  # training_gene_ids <- unique(training_gene_ids$RAP)
  # training_gene_ids <-
  #   training_gene_ids[which(training_gene_ids != "None")]
  
  universe = riceidconverter::RiceIDConvert(WEIGHTED_DF$genes, fromType = "MSU", toType = "RAP")
  universe = unique(universe$RAP)
  universe = universe[which(universe != "None")]
  
  
  results_list <- list()
  for(j in 2:ncol(WEIGHTED_DF)){
    print(j)
    #j <- 1
    if(j==2){
      candidates <-
        WEIGHTED_DF[which(WEIGHTED_DF$status ==1), ]  
    } else {
      candidates <-
        WEIGHTED_DF[which(WEIGHTED_DF$status == 0 &
                            WEIGHTED_DF[,j] == 1), ]
    }
    candidates_gene_ids <-
      riceidconverter::RiceIDConvert(candidates$genes, fromType = "MSU", toType = "RAP")
    candidates_gene_ids <- unique(candidates_gene_ids$RAP)
    candidates_gene_ids <-
      candidates_gene_ids[which(candidates_gene_ids != "None")]
    
    
 

  

  
  
  CH <- c("CANDIDATES")
  
  
  x_Hsap <-
    list(#as.character(training_gene_ids),
         as.character(candidates_gene_ids))
  names(x_Hsap) <- CH
  
  #Using as background the unique genes for the ten CH.
  
  GOterm_field <- "term_name"
  
  KEGG_LIST <- list()
  BP_LIST <- list()
  MF_LIST <- list()
  
  for (i in 1) {
    x_s <-  gprofiler2::gost(
      query = x_Hsap[[i]],
      organism = "osativa",
      ordered_query = FALSE,
      multi_query = FALSE,
      significant = TRUE,
      exclude_iea = FALSE,
      measure_underrepresentation = FALSE,
      evcodes = FALSE,
      user_threshold = 0.05,
      correction_method = "g_SCS",
      domain_scope = "annotated",
      custom_bg = as.character(universe),
      numeric_ns = "",
      sources = "KEGG",
      as_short_link = FALSE
    )
    
    
    if (!is.null(x_s$result)) {
      colnames(x_s$result)[1] <- "feature"
      x_s$result$feature <- CH[[i]]
      KEGG_LIST[[i]] <- x_s$result
      KEGG_LIST[[i]]$method
    } else {
      KEGG_LIST[[i]] <- NULL
    }
    
    
    
    
    x_s_MF <-  gprofiler2::gost(
      query = x_Hsap[[i]],
      organism = "osativa",
      ordered_query = FALSE,
      multi_query = FALSE,
      significant = TRUE,
      exclude_iea = FALSE,
      measure_underrepresentation = FALSE,
      evcodes = FALSE,
      user_threshold = 0.05,
      correction_method = "g_SCS",
      domain_scope = "annotated",
      custom_bg = as.character(universe),
      numeric_ns = "",
      sources = "GO:MF",
      as_short_link = FALSE
    )
    
    if (!is.null(x_s_MF$result)) {
      colnames(x_s_MF$result)[1] <- "feature"
      x_s_MF$result$feature <- CH[[i]]
      MF_LIST[[i]] <- reducereduntGOALL(
        df = x_s_MF$result,
        GO.ID = "term_id",
        FDR_col = "p_value",
        mode = "MF"
      )
      
    } else {
      MF_LIST[[i]] <- NULL
    }
    
    
    
    
    x_s_BP <-  gprofiler2::gost(
      query = x_Hsap,
      organism = "osativa",
      ordered_query = FALSE,
      multi_query = FALSE,
      significant = TRUE,
      exclude_iea = FALSE,
      measure_underrepresentation = FALSE,
      evcodes = FALSE,
      user_threshold = 0.05,
      correction_method = "g_SCS",
      domain_scope = "annotated",
      custom_bg = as.character(universe),
      numeric_ns = "",
      sources = "GO:BP",
      as_short_link = FALSE
    )
    
    if (!is.null(x_s_BP$result)) {
      colnames(x_s_BP$result)[1] <- "feature"
      x_s_BP$result$feature <- CH[[i]]
      BP_LIST [[i]] <- reducereduntGOALL(
        df = x_s_BP$result,
        GO.ID = "term_id",
        FDR_col = "p_value",
        mode = "BP"
      )
      
    } else {
      BP_LIST [[i]] <- NULL
    }
    
  }
  
  
  
  KEGG_LIST <- KEGG_LIST[!sapply(KEGG_LIST,is.null)]
  MF_LIST <- MF_LIST[!sapply(MF_LIST,is.null)]
  BP_LIST <- BP_LIST[!sapply(BP_LIST,is.null)]
  
  
  KEGG_LIST <- do.call(rbind,KEGG_LIST)
  MF_LIST <- do.call(rbind,MF_LIST)
  BP_LIST <- do.call(rbind,BP_LIST)
  
  
  if(!is.null(KEGG_LIST)){
    KEGG_LIST$approach <- colnames(WEIGHTED_DF)[j]
  }
  
  if(!is.null(MF_LIST)){
    MF_LIST$approach <- colnames(WEIGHTED_DF)[j]
  }
  
  if(!is.null(BP_LIST)){
    BP_LIST$approach <- colnames(WEIGHTED_DF)[j]
  }
  
  
  results <- rbind(KEGG_LIST,MF_LIST)
  results <- rbind(results,BP_LIST)
  
  
  results_list[[j-1]] <- results
  
  
  }
  message("Saving results in a XLSX file")
  

  
  results_list <- results_list[!sapply(results_list,is.null)]
  results_list <- do.call(rbind,results_list)
  
  ##names(results_list) <- colnames(WEIGHTED_DF[,-c(1)])
  openxlsx::write.xlsx(results_list, file = paste0(out_dir,"/",out_name,".xlsx"))
  message("DONE!")
  return(results_list)
  
}




  dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE/GOCOMPARE"
  out_dir <- dir
  filename = "COMPARING APPROACHES.xlsx"
  out_name <- "GO_GOCOMPARE_ALL"
  

 final_file <- ENRICHMENT_CANDIDATES_FUNC_COL(dir,filename,out_dir,out_name)
# 
# 
# dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
# out_dir <- dir
# filename = "CANDIDATES_GOCompare_ENSEMBLE_UNWEIGHTED.xlsx"
# out_name <- "GO_UNWEIGHTED"
# 
# final_file2 <- ENRICHMENT_CANDIDATES_FUNC(dir,filename,out_dir,out_name)
