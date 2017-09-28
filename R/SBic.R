## This is a generic biclustering and downstream analysis of patients by gene procedure,
## to investigate immune evasion mechanism of cancer patients
## memory.limit(size = 2000)
#library(biclust)
#args <- commandArgs(trailingOnly = TRUE)
#if (as.character(args[1]) %in% c("-h","-help","h","help","H","Help")){
  #print("Usage: Rscript Seq_Bic_Simp.R Dataset rand_seed")
  #print("Example1: Rscript Seq_Bic_Simp.R BRCA_RNAseqV2_PPITum.txt 1234")
  #print("--OR--")
  #print("Example2: Rscript Seq_Bic_Simp.R BRCA_RNAseqV2_PPITum.txt -n 1234 5")
  #quit()
#}
## Extract character string of dataset name from command line arguments
#loc_trt <- which(strsplit(as.character(args[1]), "")[[1]]==".")
#Dataname <- substr(as.character(args[1]), 1, loc_trt - 1)
#print(paste(Dataname, " is the dataset to be analyzed.", sep = ""))
#load("Data_mat.RData")
#X <- read.table(file = as.character(args[1]), header = T, row.names = 1)
#X <- as.matrix(X)
## Generic output of gene names and patient ID's
SBic <- function(X, propcutoff = 0.05, seed = 1234, outopt = "table", Dataname = "Dataset") {
  ## Next we run sequential biclustering
  ## We need to do several iterations, each one including a bicluster and a removing
  ## the bicluster with most patients
  ## Also we implement a t-test each time to the removed bicluster
  #GENE_ID<-rownames(PPI_Tum)
  if (is.matrix(X)==0 | is.numeric(X)==0){
    stop("Input data must be a numeric matrix!")
  }
  Dataset1 <- X
  Ngenes <- nrow(X)
  Npatients <- ncol(X)
  min_patient <- round(Npatients*propcutoff)
  print(paste("Minimum number of patient required for a valid cluster is ", as.character(min_patient), sep=""))
  PatientID1 <- colnames(X)
  Genename1 <- rownames(X)
  ind <- numeric()
  Data1_ColInd_rem <- Data1_RowInd_union <- numeric()
  Data1_ColInd <- Data1_RowInd <- list()
  Data1_step <- Data1_Bic <- list()
  print("Now we run sequential (bi)clustering:")
  ite <- 1
  set.seed(seed)
  ## Iteration 1:
  ## Extract the input seed
  #seed <- as.numeric(args[2])
  Data1_step[[ite]] <- biclust(Dataset1, method = BCPlaid(), shuffle = 3, row.release = 0.3, 
                               col.release = 0.3, verbose = FALSE)
  ## Obtain number of biclusters found in this step
  ncol_Bic <- numeric(Data1_step[[ite]]@Number)
  for (i in 1:Data1_step[[ite]]@Number){
    ncol_Bic[i] <- ncol(bicluster(Dataset1, Data1_step[[ite]], i)[[1]])
  }
  #This line judges whether the sequential procedure can go on!!!
  ind[ite] <- ifelse(max(ncol_Bic)<min_patient, 0, min(which(ncol_Bic>=min_patient)))
  if (ind[ite]>0){
    Data1_Bic[[ite]] <- bicluster(Dataset1, Data1_step[[ite]], ind[ite])
    Data1_ColInd[[ite]] <- which(PatientID1 %in% colnames(Data1_Bic[[ite]][[1]]))
    Data1_RowInd[[ite]] <- which(Genename1 %in% rownames(Data1_Bic[[ite]][[1]]))
    Data1_ColInd_rem <- c(Data1_ColInd_rem, Data1_ColInd[[ite]])
    #print(paste("Number of patients to remove is: ", as.character(length(Data1_ColInd_rem)),sep=""))
    Data1_RowInd_union <- c(Data1_RowInd_union,Data1_RowInd[[ite]])
    #print(paste("Number of genes involved so far is: ", as.character(length(Data1_RowInd_union)),sep=""))
    Dataset1 <- X[,-Data1_ColInd_rem]  ## Must remove from the original matrix
    print(paste("Iteration ", as.character(ite), " done!", sep=""))
    print(paste("Number of genes is ", as.character(nrow(Data1_Bic[[ite]][[1]])),
                ". Number of patients is ", as.character(ncol(Data1_Bic[[ite]][[1]])), sep = ""))
    print(paste("The remaining dataset dimension is: ", as.character(nrow(Dataset1)), " by ", 
                as.character(ncol(Dataset1)), sep = ""))
  }
  
  ## Following iterations, while loop works:
  while (ind[ite]>0){
    ite <- ite + 1
    set.seed(seed)
    Data1_step[[ite]] <- biclust(Dataset1, method = BCPlaid(), shuffle = 3, row.release = 0.3, 
                                 col.release = 0.3, verbose = FALSE)
    ncol_Bic <- numeric(Data1_step[[ite]]@Number)
    for (i in 1:Data1_step[[ite]]@Number){
      ncol_Bic[i] <- ncol(bicluster(Dataset1, Data1_step[[ite]], i)[[1]])
    }
    #This line judges whether the sequential procedure can go on!!!
    ind[ite] <- ifelse(max(ncol_Bic)<min_patient, 0, min(which(ncol_Bic>=min_patient)))
    if (ind[ite]>0){
      Data1_Bic[[ite]] <- bicluster(Dataset1, Data1_step[[ite]], ind[ite])
      Data1_ColInd[[ite]] <- which(PatientID1 %in% colnames(Data1_Bic[[ite]][[1]]))
      Data1_RowInd[[ite]] <- which(Genename1 %in% rownames(Data1_Bic[[ite]][[1]]))
      Data1_ColInd_rem <- c(Data1_ColInd_rem, Data1_ColInd[[ite]])
      Data1_RowInd_union <- c(Data1_RowInd_union, Data1_RowInd[[ite]])
      #print(paste("Number of genes involved so far is: ", as.character(length(Data1_RowInd_union)),sep=""))
      Dataset1 <- X[,-Data1_ColInd_rem]  ## Must remove from the original matrix
      print(paste("Iteration ", as.character(ite), " done!", sep = ""))
      print(paste("Number of genes is ", as.character(nrow(Data1_Bic[[ite]][[1]])),
                  ". Number of patients is ", as.character(ncol(Data1_Bic[[ite]][[1]])), sep = ""))
      print(paste("The remaining dataset dimension is: ", as.character(nrow(Dataset1)), " by ", 
                  as.character(ncol(Dataset1)), sep = ""))
    }
  }
  ## Note that the last element of ind must be 0, so we shouldn't use "num_of_Bic <- length(ind)"
  num_of_Bic <- sum(ind>0)
  print(paste("Sequential procedure finished, found ", as.character(num_of_Bic), " (bi)clusters in total!", sep = ""))
  if (outopt == "table"){
    for (i in 1:num_of_Bic){
      write.table(Data1_Bic[[i]][[1]], file = paste(Dataname, "_CL", as.character(i), ".txt", sep=""), 
                  quote = F, row.names = T, col.names = T)
    }
    print("(Bi)clusters of genes by patients output!")
  }
}
#write.table(colnames(X), file = paste(Dataname,"_PatientID",".txt",sep=""), 
            #quote = F, row.names = F, col.names = "PatientID")









