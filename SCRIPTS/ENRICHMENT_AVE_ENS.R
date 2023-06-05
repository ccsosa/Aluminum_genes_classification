  require(readxl)
  require(ggplot2)
  require(riceidconverter)
  require(gprofiler2)
  require(openxlsx)
  
  
  ENRICHMENT_CANDIDATES_FUNC_COL <- function(dir,filename,sheet,field,out_dir,out_name){
    source("https://raw.githubusercontent.com/ccsosa/reducereduntGO/master/reducereduntGO_ALL.R")
    
  
    WEIGHTED_DF <-
      as.data.frame(read_xlsx(paste0(dir, "/", filename),sheet = sheet))
    
    # training <- WEIGHTED_DF[which(WEIGHTED_DF$status == 1), ]
    # training_gene_ids <-
    #   riceidconverter::RiceIDConvert(training$genes, fromType = "MSU", toType = "RAP")
    # training_gene_ids <- unique(training_gene_ids$RAP)
    # training_gene_ids <-
    #   training_gene_ids[which(training_gene_ids != "None")]
    
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
      
   
  
    
  if(length(candidates_gene_ids)>0){
    
    
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
    
    #i <- 1
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
        #domain_scope = "annotated",
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
        #domain_scope = "annotated",
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
        #domain_scope = "annotated",
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
      KEGG_LIST$approach <- WEIGHTED_DF_datasets[[j]]#colnames(WEIGHTED_DF)[j]
    }
    
    if(!is.null(MF_LIST)){
      MF_LIST$approach <-  WEIGHTED_DF_datasets[[j]]#colnames(WEIGHTED_DF)[j]
    }
    
    if(!is.null(BP_LIST)){
      BP_LIST$approach <-  WEIGHTED_DF_datasets[[j]]#colnames(WEIGHTED_DF)[j]
    }
    
    
    results <- rbind(KEGG_LIST,MF_LIST)
    results <- rbind(results,BP_LIST)
    
    
    results_list[[j]] <- results
  } else {
    results_list[[j]] <-NULL
  }
      
      }
    message("Saving results in a XLSX file")
    
  
    
    results_list <- results_list[!sapply(results_list,is.null)]
    results_list <- do.call(rbind,results_list)
    
    ##names(results_list) <- colnames(WEIGHTED_DF[,-c(1)])
    openxlsx::write.xlsx(results_list, file = paste0(out_dir,"/",out_name,".xlsx"))
    message("DONE!")
    return(results_list)
    
  }
  
  
  
  
   #  dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE/GOCOMPARE"
   #  out_dir <- dir
   #  filename = "COMPARING APPROACHES.xlsx"
   #  out_name <- "GO_GOCOMPARE_ALL"
   #  
   # 
   # final_file <- ENRICHMENT_CANDIDATES_FUNC_COL(dir,filename,out_dir,out_name)
   # 
   
   
   # dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE/JJG"
   # out_dir <- dir
   # filename = "COMPARING APPROACHES.xlsx"
   # out_name <- "GO_GOCOMPARE_ALL"
   # 
   # 
   # final_file <- ENRICHMENT_CANDIDATES_FUNC_COL(dir,filename,out_dir,out_name)
   # 
   dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
   out_dir <- dir
   filename = "COMBINED_RESULTS.xlsx"
   sheet = "Candidate ranks"
   field = "status D1"
   out_name <- "GO_GOCOMPARE_RANKING_D1"
   
   
   final_file <- ENRICHMENT_CANDIDATES_FUNC_COL(dir,filename,sheet,field,out_dir,out_name)
   
 
 dir <- "D:/REPO_GITHUB/ALUMINUM_GENES_CLASSIFICATION/ENSEMBLE"
 out_dir <- dir
 filename = "COMBINED_RESULTS.xlsx"
 sheet = "Candidate ranks"
 field = "status D2"
 out_name <- "GO_GOCOMPARE_RANKING_D2"
 
 
 final_file <- ENRICHMENT_CANDIDATES_FUNC_COL(dir,filename,sheet,field,out_dir,out_name)
 
 