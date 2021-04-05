library(MASS)
###############################################################
# function alg is used to calculate the candidate matrices
# the function here works for both univariate and multivariate
###############################################################
alg <- function(X,Y,s=0.02,n.w,percent=1,method=c("kernel","euler"),scale="X"){  #X;Y=y;s=0.02;n.w;percent;method="euler" 
  if (scale=="Z"){
    spcovrt <- mppower(var(X),-0.5,0.000000001)  
    c.X <- scale(X,center=TRUE,scale = F)   #center X
    z <- c.X%*%spcovrt                      #standardize X
    p <- ncol(X)
    n <- nrow(X)
    q <- ncol(Y)
    
    if (method=="kernel"){   
      sig.w <- 0.1          # variance used to for kernel function as suggested in Zhu, Zeng 2006
      Mo <- matrix(0,p,p)
      if(q==1){
        SY <- (Y-mean(Y))      #center Y 
        for(i in 1:n){
          z1 <- z[i,]
          y1 <- SY[i]
          Mi <- matrix(0,p,p)
          for(j in 1:n){
            z2 <- z[j,]
            y2 <- SY[j]
            l2n <- sum(t(z1-z2)%*%(z1-z2))
            M1 <- y1*y2*exp(-sig.w*l2n/2)*sig.w*(diag(1,p)-sig.w*(z1-z2)%*%t(z1-z2))
            Mi <- Mi+M1
          }
          Mo <- Mo + Mi
        }
        CDDT.mat <- Mo/(n*n)
      } else {
        spcovrtY <- mppower(var(Y),-0.5,0.000000001)  
        c.Y <- scale(Y,center=TRUE,scale = F)   #center Y 
        #sdY <- c.Y%*%spcovrtY                     #center Y
        sdY <- c.Y
        for(i in 1:n){
          z1 <- z[i,]
          y1 <- sdY[i,]
          Mi <- matrix(0,p,p)
          for(j in 1:n){
            z2 <- z[j,]
            y2 <- sdY[j,]
            l2n <- sum(t(z1-z2)%*%(z1-z2))
            M1 <- as.vector(t(y1)%*%y2)*exp(-sig.w*l2n/2)*sig.w*(diag(1,p)-sig.w*(z1-z2)%*%t(z1-z2))
            Mi <- Mi+M1
          }
          Mo <- Mo + Mi
        }
        CDDT.mat <- Mo/(n*n)
      }
    }
    else if(method=="euler"){  
      cn.w <- n.w*percent   #actually used number of w
      sig <- sum(diag(t(X)%*%X))/n
      sig.w <- s*pi^2/sig*diag(1,p)   
      mu.w <- rep(0,p)
      matw0 <- mvrnorm(n.w,mu.w,sig.w)  # generate a n.w*p matrix, n.w w's  
      matw00 <- diag((diag(matw0%*%t(matw0))^(-0.5)), n.w, n.w)%*%matw0  # standardize w matrix
      matw <- t(matw00)
      matw <- as.matrix(matw)
      if(q==1){
        SY <- (Y-mean(Y))
        cp1 <- abs(cor(SY,cos(z%*%matw)))         #calculate the abs correlation between y and cos part
        cp2 <- abs(cor(SY,sin(z%*%matw)))
        cp <- c(cp1,cp2)                          
        cporder <- order(cp,decreasing = TRUE)     #get the index number of cp from the largest to smallest
        cporder[which(cporder > n.w)] <- cporder[which(cporder > n.w)] - n.w    
        ndx <- unique(cporder)[1:cn.w]            #get the index number of the w that satisfy the condition
        matwf <- matw[,ndx]
        C.mat <- matrix(NA,p,2* cn.w)
        for(i in 1: cn.w){
          w <- matwf[,i]
          b <- w*mean((cos(z%*%w)-mean(cos(z%*%w)))*(SY-mean(SY)))/sqrt(var(cos(z%*%w))*var(SY))     #calculate the imaginary part
          a <- -w*mean((sin(z%*%w)-mean(sin(z%*%w)))*(SY-mean(SY)))/sqrt(var(sin(z%*%w))*var(SY))   
          C.mat[,(2*i-1)] <- a
          C.mat[,(2*i)] <- b
        }
        CDDT.mat <- C.mat%*%t(C.mat)
      } else {
        spcovrtY <- mppower(var(Y),-0.5,0.000000001)  
        c.Y <- scale(Y,center=TRUE,scale = F)   #center Y 
        #sdY <- c.Y%*%spcovrtY                     #standardize Y
        sdY <- c.Y
        cp1m <- abs(cor(sdY,cos(z%*%matw)))
        cp2m <- abs(cor(sdY,sin(z%*%matw)))
        cp1 <- apply(cp1m, 2,max)         #choose the max abs correlation between y and cos part
        cp2 <- apply(cp2m, 2,max)
        cp <- c(cp1,cp2)                          
        cporder <- order(cp,decreasing = TRUE)     #get the index number of cp from the largest to smallest
        cporder[which(cporder > n.w)] <- cporder[which(cporder > n.w)] - n.w   
        ndx <- unique(cporder)[1:cn.w]            #get the index number of the w that satisfy the condition
        matwf <- matw[,ndx]
        C.mat <- matrix(0,p,2*cn.w*q)
        for(i in 1: cn.w){   #i=2
          w <- matwf[,i]
          sdcos <- matrix((cos(z%*%w)-mean(cos(z%*%w)))/as.vector(sqrt(var(cos(z%*%w)))),n,q)
          sdsin <- matrix((sin(z%*%w)-mean(sin(z%*%w)))/as.vector(sqrt(var(sin(z%*%w)))),n,q)
          a <- as.matrix(w)%*%apply(sdcos*c.Y,2,mean)    #calculate the imaginary part
          b <- -as.matrix(w)%*%apply(sdsin*c.Y,2,mean)   
          C.mat[,((2*q*i-2*q+1):(2*i*q-q))] <- a
          C.mat[,((2*i*q-q+1):(2*q*i))] <- b
        }
        CDDT.mat <- C.mat%*%t(C.mat)
      }
    } 
    else stop("Wrong argument for method=!")
  } else {
    z <- X
    p <- ncol(X)
    n <- nrow(X)
    q <- ncol(Y)
    if (method=="kernel"){   
      sig.w <- 0.1          # variance used to for kernel function as suggested in Zhu, Zeng 2006
      Mo <- matrix(0,p,p)
      if(q==1){
        SY <- (Y-mean(Y))
        #SY <- Y
        for(i in 1:n){
          z1 <- z[i,]
          y1 <- SY[i]
          Mi <- matrix(0,p,p)
          for(j in 1:n){
            z2 <- z[j,]
            y2 <- SY[j]
            l2n <- sum(t(z1-z2)%*%(z1-z2))
            M1 <- y1*y2*exp(-sig.w*l2n/2)*sig.w*(diag(1,p)-sig.w*(z1-z2)%*%t(z1-z2))
            Mi <- Mi+M1
          }
          Mo <- Mo + Mi
        }
        CDDT.mat <- Mo/(n*n)
      } else {
        spcovrtY <- mppower(var(Y),-0.5,0.000000001)  
        c.Y <- scale(Y,center=TRUE,scale = F)   #center Y 
        #sdY <- c.Y%*%spcovrtY                     #standardize Y
        sdY <- c.Y
        for(i in 1:n){
          z1 <- z[i,]
          y1 <- sdY[i,]
          Mi <- matrix(0,p,p)
          for(j in 1:n){
            z2 <- z[j,]
            y2 <- sdY[j,]
            l2n <- sum(t(z1-z2)%*%(z1-z2))
            M1 <- as.vector(t(y1)%*%y2)*exp(-sig.w*l2n/2)*sig.w*(diag(1,p)-sig.w*(z1-z2)%*%t(z1-z2))
            Mi <- Mi+M1
          }
          Mo <- Mo + Mi
        }
        CDDT.mat <- Mo/(n*n)
      }
    }
    else if(method=="euler"){  
      cn.w <- n.w*percent   #actually used number of w
      sig <- sum(diag(t(X)%*%X))/n
      #sig <- median(diag(t(X)%*%X))
      sig.w <- s*pi^2/sig*diag(1,p)
      #sig.w <- s*diag(1,p)
      mu.w <- rep(0,p)
      matw0 <- mvrnorm(n.w,mu.w,sig.w)  # generate a n.w*p matrix, n.w w's  
      matw00 <- diag((diag(matw0%*%t(matw0))^(-0.5)), n.w, n.w)%*%matw0  # normalize w matrix
      matw <- t(matw00)
      matw <- as.matrix(matw)
      if(q==1){
        #SY <- (Y-mean(Y))/as.vector(sqrt(var(Y)))
        SY <- Y-mean(Y)
        cp1 <- abs(cor(SY,cos(z%*%matw)))         #calculate the abs correlation between y and cos part
        cp2 <- abs(cor(SY,sin(z%*%matw)))
        cp <- c(cp1,cp2)                          
        cporder <- order(cp,decreasing = TRUE)     #get the index number of cp from the largest to smallest
        cporder[which(cporder > n.w)] <- cporder[which(cporder > n.w)] - n.w    
        ndx <- unique(cporder)[1:cn.w]            #get the index number of the w that satisfy the condition
        matwf <- matw[,ndx]
        C.mat <- matrix(NA,p,2* cn.w)
        for(i in 1: cn.w){
          w <- matwf[,i]
          b <- w*mean((cos(z%*%w)-mean(cos(z%*%w)))*(SY-mean(SY)))/sqrt(var(cos(z%*%w))*var(SY))     #calculate the imaginary part
          a <- w*mean((sin(z%*%w)-mean(sin(z%*%w)))*(SY-mean(SY)))/sqrt(var(sin(z%*%w))*var(SY))   
          C.mat[,(2*i-1)] <- a
          C.mat[,(2*i)] <- b
        }
        CDDT.mat <- C.mat%*%t(C.mat)
      } else {
        spcovrtY <- mppower(var(Y),-0.5,0.000000001)  
        c.Y <- scale(Y,center=TRUE,scale = F)   #center Y 
        # sdY <- c.Y%*%spcovrtY                     #standardize Y
        sdY <- scale(Y,center=TRUE,scale = F)
        cp1m <- abs(cor(sdY,cos(z%*%matw)))
        cp2m <- abs(cor(sdY,sin(z%*%matw)))
        cp1 <- apply(cp1m, 2,max)         #choose the max abs correlation between y and cos part
        cp2 <- apply(cp2m, 2,max)
        cp <- c(cp1,cp2)                          
        cporder <- order(cp,decreasing = TRUE)     #get the index number of cp from the largest to smallest
        cporder[which(cporder > n.w)] <- cporder[which(cporder > n.w)] - n.w   
        ndx <- unique(cporder)[1:cn.w]            #get the index number of the w that satisfy the condition
        matwf <- matw[,ndx]
        C.mat <- matrix(0,p,2*cn.w*q)
        for(i in 1: cn.w){   #i=2
          w <- matwf[,i]
          sdcos <- matrix((cos(z%*%w)-mean(cos(z%*%w)))/as.vector(sqrt(var(cos(z%*%w)))),n,q)
          sdsin <- matrix((sin(z%*%w)-mean(sin(z%*%w)))/as.vector(sqrt(var(sin(z%*%w)))),n,q)
          a <- as.matrix(w)%*%apply(sdcos*c.Y,2,mean)    #calculate the imaginary part
          b <- as.matrix(w)%*%apply(sdsin*c.Y,2,mean)   
          C.mat[,((2*q*i-2*q+1):(2*i*q-q))] <- a
          C.mat[,((2*i*q-q+1):(2*q*i))] <- b
        }
        CDDT.mat <- C.mat%*%t(C.mat)
      }
    } 
    else stop("Wrong argument for method=!")
  }
  
  return(CDDT.mat)
}



###########################################
# Taking power with ignoring 0 eigenvalues
###########################################
mppower = function(matrix,power,ignore){
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  m = length(eval[abs(eval)>ignore])
  tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
    t(evec[,1:m])
  return(tmp)
}
