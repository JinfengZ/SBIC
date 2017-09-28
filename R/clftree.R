#load("Data_mat.RData")
#library(biclust)
#library(rpart)
#library(rpart.plot)
## Do NOT run from here!!!
## Next we randomly sample a proportion of patients:
## Random shuffle of 1065 patients
clftree <- function(X, GroupID, CV = NULL, seed = 1234, fig = TRUE, ...){
  PatientID <- colnames(X)
  GENE_ID <- rownames(X)
  n <- ncol(X); m <- nrow(X)
  ind <- numeric()
  GroupIDcha <- character(n)
  for (i in 1:n){
    GroupIDcha[i] <- paste("IES", as.character(GroupID[i]), sep = "")
  }
  GroupIDcha[GroupIDcha=="IES0"] <- "Other"
  PPI_Tum_RP <- cbind.data.frame(t(X), GroupIDcha)
  colnames(PPI_Tum_RP)[m+1] <- "GroupID"
  Bic_RP1 <- rpart(as.factor(GroupID) ~ ., data = PPI_Tum_RP, method = "class")
  #plotcp(Bic_RP1)
  #plot(Bic_RP1,compress=F,margin=0.01)
  #text(Bic_RP1,use.n=T)
  #par(mar = c(0.01, 0.01, 0.01, 0.01))
  #prp(Bic_RP1)
  if (fig){
    png("Tree.png", width = 10, height = 8, units = "cm", res = 600)
    prp(Bic_RP1, type = 0, box.col = 5)
    dev.off()
  }
  #rpart.plot(Bic_RP1,type = 1, extra = 1, under = T, fallen.leaves = F, digits = 2, 
             #varlen = 0, faclen = 0, cex = NULL, tweak = 1.4, snip = FALSE, box.palette = 0, 
             #shadow.col = 0, col = "darkblue")
  #printcp(Bic_RP1)
  over_err <- sum(predict(Bic_RP1, type = "class")!=PPI_Tum_RP[,m+1])/n #Overall error rate
  print(paste("Overall prediction error rate is ", as.character(over_err), ".", sep = ""))
  Class_pred <- predict(Bic_RP1, type = "class"); Class_true <- PPI_Tum_RP[,m+1]
  print("Confusion matrix of prediction is as follows:")
  print(table(Class_pred, Class_true))
  #Number of folds, subject to change
  if (is.null(CV)==1){
    print("User decide not to do cross validation, although that's highly recommended.")
    return(list(Bic_RP1, rep(0,5)))
    quit()
  }
  set.seed(seed)
  size_per_fold <- n%/%CV
  # Use aaa as the randomization scheme
  # block is the randomized subgroup index. Each observation is assigned a value between {1,2,3,...,K}
  aaa <- rank(runif(n))
  block <- aaa%%CV + 1
  #block<-as.factor(block)
  Err_clf <- numeric(CV)
  Mod_clf <- Pred_clf <- list(CV)
  #This algorithm implements K-fold CV and computes average error rates
  for (k in 1:CV){
    #data_to_use<-XY[block!=k,]
    Mod_clf[[k]] <- rpart(as.factor(GroupID) ~ ., data = PPI_Tum_RP[block!=k,], method = "class")
    Pred_clf[[k]] <- predict(Mod_clf[[k]], newdata = PPI_Tum_RP[block==k,], type = "class")
    Err_clf[k] <- sum(Pred_clf[[k]]!=PPI_Tum_RP[block==k,m+1])/sum(block==k)
  }
  return(list(Bic_RP1, Err_clf))
}


