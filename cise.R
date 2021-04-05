###############################################################################
#function to get the CISE sparse estimation
#the code here is rewritten based on the matlab code by Chen, Zou and Cook(2012)
#we also modify it to work for our method
#extend the original code to work for save and phd method as well
#package pracma has to be installed
###############################################################################
library(pracma)
source("algorithm.R")
cise <- function(y,x,di,pw,method){ #
  n <- dim(x)[1]
  p <- dim(x)[2]
  d <- di
  nslice <- 6
  if (method=="SIR"){
    sir <- SIR(y,x,nslice)
    M <- sir$M
    N <- sir$N
  } else if(method=="PFC"){
    pfc <- PFC(y,x)
    M <- pfc$M
    N <- pfc$N
  }else if(method=="EULER"){
    euler <- EULER(X=x,Y=y,n.w=1000,percent=0.2)
    M <- euler$M
    N <- euler$N
  }else if(method=="KERNEL"){
    kernel <- KERNEL(X=x,Y=y)
    M <- kernel$M
    N <- kernel$N
  }else if(method=="SAVE"){
    save <- SAVE(x,y,nslice)
    M <- save$M
    N <- save$N
  }else if(method=="PHD"){
    phd <- PHD(x,y)
    M <- phd$M
    N <- phd$N
  }
  else stop("Wrong argument of method=!")
  
  nin <- length(pw)
  stop <- c()
  for(i in 1:nin){ #i=91
    f <- biccis(pw[i],M,N,d,p,n)$fv
    stop[i] <- f
  }
  ind <- which.min(stop)
  lap <- ind
  dlap <- lap
  bicresult <- biccis(dlap,M,N,d,p,n)
  return(list(fv=bicresult$fv,beta=bicresult$beta,st=bicresult$st))
}
# This is the main function that achives coordinate-independent sparse estimation using
# different method to calculat M and N matrices. 
# input
#  y: response vector.
#  x: predictors matrix.
#  di: the dimension of the central subspace
#  pw: the range of the penalty parameter of \lambda is from 0 to pw. 
#  method: model to calculate M and N matrices. 
# output
#  fv: the value of the objective function at the estimator
#  beta:  the estimator of the central subspace based on BIC criterion.
#  st: a vector with either 0 or 1 element to indicate which variable is to be selected by CISE.


SIR <- function(y,x,nslice){  #checked
  n <- dim(x)[1]
  p <- dim(x)[2]
  xo <- scale(x,center = TRUE, scale = FALSE)
  N <- t(xo)%*%xo/n
  index <- order(y)
  cmean <-apply(x,2,mean)
  group <- sort((1:n)%%nslice)
  sigmaeta <- matrix(0,p,p)
  diffmean <- rep(1,p)
  for(i in 0:(nslice-1)){
    nsindx <- index[group==i]
    diffmean <- colMeans(x[nsindx,]) - cmean
    sigmaeta = sigmaeta + length(nsindx)/n*diffmean%*%t(diffmean)
  }
  M = (sigmaeta + t(sigmaeta))/2
  return(list(M=M,N=N))
}

SAVE = function(x,y,nslice){
  spcovrt <- mppower(var(x),-0.5,0.000000001)  #percent is the percentile want to keep
  c.X <- scale(x,center=TRUE,scale = F)   #center X
  z <- c.X%*%spcovrt 
  n = nrow(x)
  p = ncol(x)
  H=nslice
  N <- t(c.X)%*%c.X/n
  sliced_y = numeric(n)
  slice_prop = numeric(H)
  divpt = c(min(y)-0.1,min(y)+((range(y)[2]-range(y)[1])/H)*c(1:H))
  save_matrix = matrix(0,p,p)
  for(i in 1:H){
    sliced_y[(y<=divpt[i+1])&(y>divpt[i])]=i
    slice_prop[i] = sum(sliced_y[sliced_y==i])/(n*i)
    slice_i = matrix(z[sliced_y==i,],length(z[sliced_y==i,])/p,p)
    if (slice_prop[i]!=0)
      save_matrix = save_matrix + slice_prop[i]*mppower((diag(1,p)- t(slice_i)%*%slice_i/(length(slice_i)/p)),2,0.000000001)
    
  }
  varXhf <- mppower(var(x),0.5,0.000000001)
  M <- varXhf%*%save_matrix%*%varXhf
  return(list(M=save_matrix,N=N))
}

PHD = function(x,y){
  r = residuals(lm(y~x))
  xbar=colMeans(x)
  n = nrow(x)
  p = ncol(x)
  xo <- scale(x,center = TRUE, scale = FALSE)
  N <- t(xo)%*%xo/n
  sigma_rxx = matrix(0,p,p)
  for(i in 1:n){
    sigma_rxx = sigma_rxx +r[i]*t(t(x[i,]-xbar))%*%t(x[i,]-xbar)
  }
  m_phd = solve(cov(x))%*%(sigma_rxx/n)
  M <- m_phd%*%t(m_phd)
  return(list(M=M,N=N))
}

PFC <- function(y,x){  #checked
  n <- dim(x)[1]
  p <- dim(x)[2]
  xo <- scale(x,center = TRUE, scale = FALSE)
  N <- t(xo)%*%xo/n
  
  yy <- cbind(sqrt(abs(y)), y, y^2)
  yf <- scale(yy,center = TRUE,scale = FALSE)
  Pf <- yf%*%solve(t(yf)%*%yf)%*%t(yf)
  Xf <- Pf%*%xo
  Sigmaf <- t(Xf)%*%Xf/n
  
  M <- (Sigmaf + t(Sigmaf))/2
  return(list(M=M,N=N))
} 


KERNEL <- function(X,Y){
  CDDT.mat <- alg(X,Y,n.w=1000,method="kernel")
  c.X <- scale(X,center=TRUE,scale = F)   
  N <- t(c.X)%*%c.X/nrow(X)
  return(list(M=CDDT.mat,N=N))
}



EULER <- function(X,Y,s=0.02,n.w=1000,percent=0.2){
  CDDT.mat <- alg(X,Y,n.w=1000,method="euler")
  c.X <- scale(X,center=TRUE,scale = F)   
  N <- t(c.X)%*%c.X/nrow(X)
  return(list(M=CDDT.mat,N=N))
}



biccis <- function(lam,M,N,d,p,n){ #lam=i/100-0.01;lam=-0.009
  M0 <- M
  N0 <- N
  if(lam < 0){
    f=Inf
  } else{
    lambda = lam #lam=0
    mysqrtr <- mysqrt(N)
    #A <- mysqrtr$A
    B <- mysqrtr$B
    G <- B%*%M%*%B
    egG <- eigen(G)
    a2 <- egG$vectors
    b2 <- egG$values
    b1 <- order(b2)
    if(d==2){          ##only consider d=1,2 based on the original code
      xo <- B%*%cbind(a2[,b1[p]],a2[,b1[p-1]])   
    } else if(d==1){
      xo <- B%*%a2[,b1[p]]
    }
    ecisr <- ecis(xo,lambda,p,M,N)
    beta <- ecisr$beta
    st <- ecisr$st
    M <- M0
    N <- N0
    f <- -sum(diag(t(beta)%*%M%*%beta)) + d*(sum(st)-d)*log(n)/n
  }
  return(list(fv=f,beta=beta,st=st))
} 


mysqrt <- function(N){  
  eM <- eigen(N)
  a1 <- eM$vectors
  b1 <- eM$values
  if(length(a1)==1){
    x <- a1*sqrt(b1)%*%a1
    y <- a1*1/sqrt(b1)*a1
  } else{
    x <- a1%*%diag(sqrt(b1))%*%t(a1)
    x <- (x+t(x))/2
    y <- a1%*%diag(1/sqrt(b1))%*%t(a1)
    y <- (y+t(y))/2
  }
  return(list(A=x,B=y))
}  

mywt <- function(b){  
  p0 <- max(dim(b))
  x <- c()
  for (i in 1:p0){
    x[i] <- sqrt(sum(b[i,]^2))^(-0.5)  #norm take power -0.5?? check the paper
  }
  return(x=x)
}


ecis <- function(xo,lambda,p,M,N){ #x=xo, 
  if(lambda==0){
    beta=xo
    st=rep(1,p)
  } else{
    beta=matrix(0,dim(xo)[1],dim(xo)[2])
    d <- dim(xo)[2]
    b0 <- xo
    the <- mywt(xo)
    dis <- 1
    st <- rep(1,p)
    H0 <- diag(the)%*%diag(1/sqrt(diag(b0%*%t(b0))))
    
    B <- mysqrt(N)$B
    G <- B%*%M%*%B
    k=0
    while (dis > 10^(-6) ) {
      k=k+1
      weig <- eigen(G-0.5*lambda*B%*%H0%*%B) #dim(H0)
      td <- weig$vectors
      ty <- weig$values
      tp <- max(dim(H0))
      b1 <- order(ty)
      if(d==2){
        b=B%*%cbind(td[,b1[tp]],td[,b1[tp-1]])
      } else if(d==1){
        b=B%*%td[,b1[tp]]
      }
      dis <- subspace(b0,b)
      st0 <- st
      delpr <- delp(b,st,p,M,N) 
      N <- delpr$N
      M <- delpr$M
      b <- delpr$x
      st <- delpr$y
      ix <- delpr$ix
      if(sum(st) < sum(st0)){
        the <- the[-ix] 
        A <- mysqrt(N)$A
        B <- mysqrt(N)$B
        G <- B%*%M%*%B
      }
      if(length(the)==1){
        H0 <- the*1/abs(b)
        b0 <- b
        break
      }else{
        H0 <- diag(the)%*%diag(1/sqrt(diag(b%*%t(b))))   
        b0 <- b
        if(dim(H0)[1] <=2){
          break
        }
      }
    }
    te <- 0
    for (i in 1:p) {
      if(st[i]==1){
        te <- te+1
        b <- as.matrix(b)
        beta[i,(1:d)] <- b[te,(1:d)] #dim(beta)
      }else if(st[i]==0){beta[i,(1:d)] <- rep(0,d)}
    }
    
  }
  return(list(beta=beta,st=st))
}


delp <- function(b,st,p,M,N){
  t=0
  s=0
  y=st
  #b <- round(b,5)
  ix <- c()
  for(i in 1:p){
    if(st[i]!=0){
      t=t+1
      if(sqrt(sum(b[t,]^2)) <= 10^(-6)){
        s=s+1
        y[i]=0
        ix[s]=t
      }
    }
  }
  if(sum(y) < sum(st)){
    b <- b[-ix,]
    M <- M[-ix,-ix]
    N <- N[-ix,-ix]
    x = b
  } else {
    ix=0
    x=b
  }
  return(list(x=x,y=y,ix=ix,N=N,M=M))
}


