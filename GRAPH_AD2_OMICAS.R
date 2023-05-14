require(igraph);require(riceidconverter);require(DESeq2);require(WGCNA);require(parallel);#equire(HiClimR)
dir <- "/users/ccsosaa/pecanpy"
setwd(dir)

setwd(dir)

#usethis::edit_r_environ()
graph_base <- as.matrix(read.table(paste0(dir,"/","RICENETPPI.tsv"),header = T))
x2 <- igraph::graph_from_edgelist(graph_base,directed = FALSE)
components <- igraph::clusters(x2, mode="weak")
biggest_cluster_id <- which.max(components$csize)
# ids
vert_ids <- V(x2)[components$membership == biggest_cluster_id]
# subgraph
x2 <- igraph::induced_subgraph(x2, vert_ids)
vertices <- as.data.frame(as.matrix(V(x2)))
x_rice <- riceidconverter::RiceIDConvert(myID = row.names(vertices),fromType = "MSU",toType = "RAP")
vertices$name <- row.names(vertices);
vertices$RAP <- NA
vertices$counts <- NA
vertices$CHECK <- NA


#i <- 10197
for(i in 1:nrow(vertices)){
  # print(i)
  #i <- 174#0197
  x_i <- x_rice[which(x_rice$MSU==vertices$name [[i]]),2]
  x_i <- x_i[which(x_i!="None")]
  x_i <- as.character(x_i[which(!is.na(x_i))])
  if(length(x_i)>1){
    vertices$CHECK[[i]] <- "CHECK (SEVERAL GENES)"
    vertices$RAP[[i]] <- paste(x_i,collapse = "//")
    vertices$counts[[i]] <- length(x_i)
  } else if(length(x_i)==1){
    vertices$CHECK[[i]] <- "OK (ONE GENE)"
    vertices$RAP[[i]] <- x_i
    vertices$counts[[i]] <- 1
  } else if(length(x_i)==0){
    vertices$CHECK[[i]] <- "OK (NO GENE)"
    vertices$RAP[[i]] <- NA
    vertices$counts[[i]] <- 0
  }
  
};rm(i)
write.table(vertices,paste0(dir,"/","Vertices.tsv"),na = "",row.names = F,quote = F,sep = "\t")
#counts in percentage
round((tapply(vertices$counts,vertices$counts,length)/17228)*100,2)
#counts in number
tapply(vertices$counts,vertices$counts,length)

V(x2)$RAP <- vertices$RAP
V(x2)$counts <- vertices$counts
V(x2)$CHECK <- vertices$CHECK

igraph::write.graph(x2,file = paste0(dir,"/","BIG_COMP.graphml"),format = "graphml")
igraph::write.graph(x2,file = paste0(dir,"/","BIG_COMP.edge"),format = "ncol")
x2_a <- read.table(paste0(dir,"/","BIG_COMP.edge"),header = F)
write.table(x2_a,paste0(dir,"/","BIG_COMP_s.edge"),na = "",quote = F,col.names = F,row.names = F,sep = "\t")

edges <- x2_a
edges$weight <- NA
################################################################################
#LOAD WGCNA
atlas_file <- read.table("/users/ccsosaa/pecanpy/JJG_DEG/COUNTS/COUNTS.tsv",header = T,row.names = 1)
genes <- row.names(atlas_file)
atlas_file <- t(atlas_file);colnames(atlas_file) <- genes

datExpr0 <- as.data.frame(log2(atlas_file+1))

#require(WGCNA)

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
################################################################################
cl <- parallel::makeCluster(18)
parallel::clusterExport(cl, varlist=c("datExpr0","atlas_file","edges"),envir=environment())


edges2 <- edges
edges2$weight <- NA
system.time(
  
  # # system.time(
  weight <- parallel::parLapplyLB(cl,
                                  X = seq_len(nrow(edges)),
                                  #X=seq_len(5000),
                                  fun = function (i){ 
                                    #for(i in 1:2500){#nrow(edges)){
                                    i <- 1
                                    print(i)
                                    #utils::setTxtProgressBar(pb, i)
                                    
                                    x_1 <- riceidconverter::RiceIDConvert(edges$V1[i],fromType = "MSU",toType = "RAP")
                                    x_2 <- riceidconverter::RiceIDConvert(edges$V2[i],fromType = "MSU",toType = "RAP")
                                    
                                    x_1_row <- as.character(x_1$RAP) %in% colnames(atlas_file)#x_Cor_VEC$row
                                    x_2_row <- as.character(x_2$RAP) %in% colnames(atlas_file)#x_Cor_VEC$col
                                    
                                    x_1 <- x_1$RAP[as.character(x_1$RAP) %in% colnames(atlas_file)]#x_Cor_VEC$row]
                                    x_2 <- x_2$RAP[as.character(x_2$RAP) %in%  colnames(atlas_file)]#x_Cor_VEC$col]
                                    
                                    x_1_row <- as.character(x_1) %in%  colnames(atlas_file)#x_Cor_VEC$row
                                    x_2_row <- as.character(x_2) %in%  colnames(atlas_file)#x_Cor_VEC$col
                                    
                                    
                                    if(length(x_1_row)==0){
                                      x_1_row <- FALSE
                                    }
                                    
                                    if(length(x_2_row)==0){
                                      x_2_row <- FALSE
                                    }
                                    
                                    #edges2$weight <- weight  
                                    if(length(x_1_row)>1){
                                      x_cor_1 <- as.matrix(datExpr0[,colnames(datExpr0) %in% as.character(x_1)])
                                      if(ncol(x_cor_1)>1){
                                        x_cor_1 <- as.matrix(rowMeans(datExpr0[,colnames(datExpr0) %in% as.character(x_1)],na.rm = T))
                                      }
                                    } else {
                                      x_cor_1 <- as.matrix(datExpr0[,colnames(datExpr0) %in% as.character(x_1)])
                                    }
                                    
                                    if(length(x_2_row)>1){
                                      x_cor_2 <- as.matrix(datExpr0[,colnames(datExpr0) %in% as.character(x_2)])
                                      if(ncol(x_cor_2)>1){
                                        x_cor_2 <- as.matrix(rowMeans(datExpr0[,colnames(datExpr0) %in% as.character(x_2)],na.rm = T))
                                      }
                                    } else {
                                      x_cor_2 <- as.matrix(datExpr0[,colnames(datExpr0) %in% as.character(x_2)])
                                    }
                                    
                                    if(sum(x_cor_1,na.rm = T)==0){
                                      x_cor_1 <- NULL
                                    } 
                                    
                                    if(sum(x_cor_2,na.rm = T)==0){
                                      x_cor_2 <- NULL
                                    } 
                                    
                                    if(is.null(x_cor_1) | is.null(x_cor_2)){
                                      xW <- 0
                                    } else if(nrow(x_cor_1)>0 & nrow(x_cor_2)>0){
                                      xW <- stats::cor(cbind(x_cor_1,x_cor_2))[2,1]
                                      #xW <- HiClimR::fastCor(cbind(x_cor_1,x_cor_2),upperTri = T,verbose = F,optBLAS = T)[2,1]
                                    } else {
                                      xW <- 0
                                    }
                                    #rm(x_cor_1,x_cor_2,x_1,x_2,x_1_row,x_2_row)
                                    xW
                                  })
)
parallel::stopCluster(cl)


edges2$weight <- as.numeric(unlist(weight))

write.table(edges2,paste0(dir,"/","BIG_COMP_W.edge"),na = "",row.names = F,col.names = F,quote = F,sep = "\t")
