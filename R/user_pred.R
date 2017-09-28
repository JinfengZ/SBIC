## This is a handy script that takes the gene expressions as input and outputs the patient cluster
## based on a fitted tree model
## XXX.txt can be a new data set with only the classifier genes but same or different patients,
## or even a simulated data, but it needs to have row (gene) names and column (patient) ID's
#args <- commandArgs(trailingOnly = TRUE)
#if (as.character(args[1]) %in% c("-h","-help","h","help","H","Help")){
  #print("Usage: Rscript User_pred_Script.R Dataset")
  #print("Example1: Rscript User_pred_Script.R XXX.txt")
  #quit()
#}
#load("Immuno_table.RData")
#loc_ <- which(strsplit(as.character(args[1]),"")[[1]]==".")
#Dataname_ <- substr(as.character(args[1]), 1, loc_ - 1)
#user_data <- t(read.table(file = as.character(args[1]), header = T, row.names = 1))
#user_data<-t(read.table(file="data_sim.txt",header=T,row.names=1))
## tree_overall is the R object of the tree fitted by rpart.
## Of course you need to have it saved in advance into a work space. The object name is "Bic_RP"
#load("Tree.RData")
#clf_gene <- c("IL2RG","ABCB1","DCN","CD40L","LCK","SELP","PD1","G6PD","ESR1","CD3G","BAX","CCR5")
#cut_off <- c(9.3,7,2.7,14,7.6,4.6,10,11,3.8,10,2.2,7.7,7.1)
#if (length( intersect(colnames(user_data),clf_gene) ) < 12){
  #print("Not all classifier genes are present. Please check data!")
  #quit()
#}

user_pred <- function(X, object, mechanism = immuno_table, Dataname = "data_sim", ...){
  #gene_ind <- find_ind(clf_gene,gene_name_simp)[[2]]
  #user_pred<-predict(Bic_RP, newdata=as.data.frame(t(user_data)), type="class")
  #full_mat <- matrix(data = NA, nrow = nrow(X), ncol = ncol(X))
  #full_mat[,gene_ind] <- user_data
  #rownames(full_mat) <- paste("Patient", as.character(1:1065), sep = "")
  #colnames(full_mat) <- gene_name_simp
  #full_mat <- as.data.frame(full_mat)
  mecha_eva <- imm_mech <- character()
  #find_ind(colnames(user_data), clf_gene_input[1:12])[[2]]
  user_pred1 <- predict(object = object, newdata = as.data.frame(t(X)), type = "class")
  user_pred1 <- as.character(user_pred1)
  #table(user_pred)
  for (i in 1:length(user_pred1)){
    mecha_eva[i] <- as.character(mechanism[which(mechanism[,1]==user_pred1[i]),2])
    imm_mech[i] <- as.character(mechanism[which(mechanism[,1]==user_pred1[i]),3])
  }
  ## Output results:
  user_output <- cbind.data.frame(colnames(X), user_pred1, mecha_eva, imm_mech)
  colnames(user_output) <- c("Patient ID", "IES", "Mechanism", "Suggested Immunotherapies")
  write.table(user_output, file = paste(Dataname, "_pred.csv", sep = ""), quote = F, row.names = F, sep = ",")
  return(user_output)
}
