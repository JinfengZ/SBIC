ttest_general <- function (X, classlabel) 
{
  ngenes <- nrow(X)
  CL <- as.numeric(as.factor(classlabel)) - 1
  if (is.numeric(X) == 0 | (is.matrix(X) == 0 & is.data.frame(X) == 0)) 
    stop("Error: Input data must be a numeric matrix or data frame!")
  if (length(unique(CL)) != 2) 
    stop("Error: Input data must have 2 phenotypic groups!")
  Xa <- X[, CL == 0]
  Xb <- X[, CL == 1]
  teststat <- rawp <- nmiss <- numeric(ngenes)
  for (i in 1:ngenes) {
    myttest <- t.test(Xa[i, ], Xb[i, ], alternative = "two.sided", 
                      mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)
    teststat[i] <- myttest$statistic
    rawp[i] <- myttest$p.value
    nmiss[i] <- sum(is.na(X[i,]))
  }
  rr <- data.frame(teststat = teststat, rawp = rawp)
  colnames(rr) <- c("teststat", "rawp")
  for (i in 1:ngenes){
    if (nmiss[i] > 0) {
      print(paste("Gene No.", as.character(i), " contains missing values.", " Please check data!", sep = ""))
    }
  }
  if (sum(nmiss) > 0) {
    print("Differential expression analysis implemented, ignoring the missing values above!")
  }
  return(rr)
}


#X.test <- matrix(data = rnorm(10000), nrow = 100, ncol = 100)
#rownames(X.test) <- paste("Gene", as.character(1:100), sep = "")
#colnames(X.test) <- paste("Sample", as.character(1:100), sep = "")
#for (i in 1:20){
  ## First 20 genes are DE, signal strengths somewhere between 0.5 and 1.5
  #X.test[i, 1:50] <- X.test[i, 1:50] + runif(1, 0.5, 1.5)
#}
#lab.test1 <- c(rep(1,50), rep(-1,50))
#lab.test2 <- c(rep(-1,50), rep(1,50))
#lab.test3 <- c(rep("puppy",50), rep("cat",50))
#lab.test4 <- c(rep("normal",50), rep("tumor",50))
#my.test1 <- ttest_general(X = X.test, classlabel = lab.test1)
#my.test2 <- t_test_genes(X = X.test, classlabel = lab.test2)
#my.test3 <- t_test_genes(X = X.test, classlabel = lab.test3)
#my.test4 <- t_test_genes(X = X.test, classlabel = lab.test4)
## Think about why my.test2$teststat==my.test4$teststat and my.test1$teststat==my.test3$teststat
## Hint: Please try as.numeric(as.factor(...)) - 1 and see why:)

#X.test[1,1] <- X.test[3,10] <- X.test[56,78] <- X.test[100,4] <- X.test[100,69] <- NA
#tt <- rnorm(20, 0, 1)
#tt[1] <- NA
#myttest <- t.test(tt[1:10], tt[11:20], alternative = "two.sided", 
                  #mu = 0, paired = FALSE, var.equal = FALSE, conf.level = 0.95)