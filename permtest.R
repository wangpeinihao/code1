########################################################
#function used to determine the structure dimension
#the code here is a modification of Wang and Yin(2015)
#the original code is for matlab
#######################################################
source("algorithm.R")
permtest <- function(X,Y,n.w=1000,percent=0.2,way=c("kernel","euler"), alpha=0.05,ptime){
  n <- dim(X)[1]
  p <- dim(X)[2]
  Mhat <- alg(X=X,Y=Y,n.w=n.w,percent=percent,method=way,scale = "Z")
  dcom <- svd(Mhat)
  dU <- dcom$u
  dd <- dcom$d
  m <- 0
  pvalue=0
  while ((pvalue < alpha) && (m < p-1) ) {
    Lamdahatm <- dd[m+1]-1/(p-m-1)*sum(dd[(m+2):p]) #test statistic
    Lamdaperm <- c(0,ptime)                                 #to store the test statistic in for Loop
    for(pi in 1:ptime){ 
      Xnew = X%*%dU
      Ynew=Y
      Xdata = Xnew
      or=sample(n,n)
      Xnew[,1:m]=Xdata[,1:m]
      Xnew[,(m+1):p]=Xdata[or,(m+1):p]
      Mhatpi <- alg(X=Xnew,Y=Ynew,n.w=n.w,percent=percent,method=way)
      permd <- svd(Mhatpi)$d
      Lamdaperm[pi] <- permd[m+1]-1/(p-m-1)*sum(permd[(m+2):p])
    }
    pvalue <- sum(Lamdaperm > Lamdahatm)
    m=m+1
  }
  estd <- m-1
  return(estd)
}

