###############################################################################
#function to get the Sequential sparse estimation
#the code here is revised based on Yin and Hilafu's code, Yin and Hilafu(2014)
#we also modify it to work for our method
#Yin and Hilafu's code is a modification of Li's code(2007) and use Hui Zou's partial code
#package pracma/lars/Rlab has to be installed
###############################################################################
source("permtest.R")
source("algorithm.R")
center<-function(v)  v - mean(v)
norm<-function(v)  
{ 
  sumv2<-sum(v^2)
  if(sumv2 == 0) sumv2<-1
  v/sqrt(sumv2)
}
orthnormal<-function(X)
{
  X<-as.matrix(X)
  n<-nrow(X)
  p<-ncol(X)
  W<-NULL
  if(p > 1) {
    W<-cbind(W, X[,1])
    for(k in 2:p) {
      gw<-rep(0, n)
      for(i in 1:(k-1)) {
        gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
        gw<-gw + gki * W[,i]
      }
      W<-cbind(W, X[,k] - gw)
    }
  } else {
    W<-cbind(W, X[,1])
  }
  
  W<-apply(W, 2, norm)
  W
}
cov.x<-function(X)
{
  Xc<-apply(X, 2, center)
  t(Xc) %*% Xc / nrow(Xc)
}

mat.sqrt<-function(A)
{
  ei<-eigen(A)
  d<-ei$values
  d<-(d+abs(d))/2
  d2<-sqrt(d)
  ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(ans)
}

mat.sqrt.inv<-function(A)
{
  ei<-eigen(A)
  d<-ei$values
  d<-(d+abs(d))/2
  d2<-1 / sqrt(d)
  d2[d == 0]<-0
  ans<-ei$vectors %*% diag(d2) %*% t(ei$vectors)
  return(ans)
}
######################################################################################
# functions for sparse principal components wrt I inner product 				 #
# (copied from Hui Zou's elasticnet library)							 #
######################################################################################

updateRR<-
  function(xnew, R = NULL, xold, lambda,eps = .Machine$double.eps)
  {
    xtx <- (sum(xnew^2)+lambda)/(1+lambda)
    norm.xnew <- sqrt(xtx)
    if(is.null(R)) {
      R <- matrix(norm.xnew, 1, 1)
      attr(R, "rank") <- 1
      return(R)
    }
    Xtx <- drop(t(xnew) %*% xold)/(1+lambda)
    r <- backsolvet(R, Xtx)
    rpp <- norm.xnew^2 - sum(r^2)
    rank <- attr(R, "rank")	
    if(rpp <= eps)
      rpp <- eps
    else {
      rpp <- sqrt(rpp)
      rank <- rank + 1
    }
    R <- cbind(rbind(R, 0), c(r, rpp))
    attr(R, "rank") <- rank
    R
  }

solvebeta<-function(x, y, paras, max.steps, sparse=c("penalty","varnum"), eps = .Machine$double.eps)
{
  sparse <- match.arg(sparse)
  if (missing(sparse)) sparse <- "penalty"
  nm <- dim(x)
  n <- nm[1]
  m <- nm[2]
  im <- seq(m)
  one <- rep(1, n)
  vn <- dimnames(x)[[2]]
  lambda<-paras[1]
  if(lambda>0){
    maxvars <- m
  }
  if (lambda==0) {
    maxvars <- min(m,n-1)
    if (m==n){
      maxvars<-m
    }
  }
  d1 <- sqrt(lambda)
  d2 <- 1/sqrt(1 + lambda)
  Cvec <- drop(t(y) %*% x) * d2
  ssy <- sum(y^2)
  residuals <- c(y, rep(0, m))
  if(missing(max.steps)) {max.steps <- 50 * maxvars}
  penalty<-max(abs(Cvec))
  if (sparse=="penalty" && penalty*2/d2<=paras[2])
  { 
    beta<-rep(0,m)
  }
  else {
    beta <- rep(0, m)
    first.in <- integer(m)
    active <- NULL
    ignores <- NULL
    actions <- as.list(seq(max.steps))
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    while((k < max.steps) & (length(active) < maxvars - length(ignores))) 
    {
      action <- NULL
      k <- k + 1
      inactive <- if(k == 1) im else im[ - c(active, ignores)]
      C <- Cvec[inactive]
      Cmax <- max(abs(C))
      if(!any(drops)) {
        new <- abs(C) == Cmax 
        C <- C[!new]
        new <- inactive[new]       
        for(inew in new) {
          R <- updateRR(x[, inew], R, x[, active], lambda) 
          if(attr(R, "rank") == length(active)) {
            nR <- seq(length(active))
            R <- R[nR, nR, drop = FALSE]
            attr(R, "rank") <- length(active)
            ignores <- c(ignores, inew)
            action <- c(action,  - inew)
          }
          else {
            if(first.in[inew] == 0)
              first.in[inew] <- k
            active <- c(active, inew)
            Sign <- c(Sign, sign(Cvec[inew]))
            action <- c(action, inew)
          }
        }
      }
      else action <-  - dropid
      Gi1 <- backsolve(R, backsolvet(R, Sign))
      A <- 1/sqrt(sum(Gi1 * Sign))
      w <- A * Gi1
      u1<-drop(x[,active,drop=FALSE]%*%w*d2)  
      u2<-rep(0,m)
      u2[active]<-d1*d2*w
      u<-c(u1,u2)
      if(lambda>0){
        maxvars <- m-length(ignores)
      }
      if (lambda==0){
        maxvars <- min(m-length(ignores),n-1)
      }
      if(length(active) == maxvars - length(ignores)) {
        gamhat <- Cmax/A
      }
      else {
        a <- (drop(u1 %*% x[,  - c(active, ignores)]) + d1 *
                u2[ - c(active, ignores)]) * d2
        gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
        gamhat <- min(gam[gam > eps], Cmax/A)                        
        Cdrop <- c(C - gamhat * a,  - C + gamhat * a) - (Cmax - gamhat * A)
      }
      dropid <- NULL
      b1 <- beta[active]
      z1 <-  - b1/w
      zmin <- min(z1[z1 > eps], gamhat)
      if(zmin < gamhat) {
        gamhat <- zmin
        drops <- z1 == zmin
      }
      else drops <- FALSE
      beta2<-beta          
      beta[active] <- beta[active] + gamhat * w
      residuals <- residuals - (gamhat*u)
      Cvec <- (drop(t(residuals[1:n]) %*% x) + d1 * residuals[ - (1:n)]) * d2
      penalty <- c(penalty,penalty[k]-abs(gamhat*A))
      if(sparse=="penalty" && rev(penalty)[1]*2/d2<=paras[2]){
        s1<-rev(penalty)[1]*2/d2
        s2<-rev(penalty)[2]*2/d2
        beta<-(s2-paras[2])/(s2-s1)*beta+(paras[2]-s1)/(s2-s1)*beta2
        beta<-beta*d2
        break
      }
      if(any(drops)) {
        dropid <- seq(drops)[drops]
        for(id in rev(dropid)) {
          R <- downdateR(R, id)
        }
        dropid <- active[drops]
        beta[dropid] <- 0
        active <- active[!drops]
        Sign <- Sign[!drops]
      }
      if(sparse=="varnum" && length(active)>=paras[2]){
        break
      }
      if(!is.null(vn))
        names(action) <- vn[abs(action)]
      actions[[k]] <- action 
    }
  }
  
  return (beta)
}


###################
#for univariate sir
###################
comp.sdr.sir<-function(X, y, h, d.true)
{
  n = nrow(X)
  p = ncol(X)
  Sigma.x = cov.x(X)
  Sigma.x2 = mat.sqrt(Sigma.x)
  Sigma.x.inv2 = mat.sqrt.inv(Sigma.x)
  
  # standardize X
  Z = apply(X, 2, center) %*% Sigma.x.inv2
  
  data = cbind(Z,y)
  data = data[order(data[,p+1]), ]
  
  znew = data[,1:ncol(X)]
  nx = nrow(X)/h
  zh = matrix(0, ncol(X),h)
  for (sl in 1:h)  
  {
    znewh = znew[((sl-1)*nx+1):(sl*nx), ]
    zh[ , sl] = colMeans(znewh)
  }
  M.sir.z = zh%*%t(zh)/h
  M.sir = Sigma.x2 %*% M.sir.z %*% Sigma.x2
  
  d = d.true
  
  m = mat.sqrt(M.sir)
  v = eigen(M.sir.z)$vectors[,1:d]
  if(d == 1) v = matrix(v, ncol=1)
  beta.sdr = Sigma.x.inv2 %*% v
  beta.sdr = apply(beta.sdr, 2, norm)
  ans = list(m=m, G=Sigma.x, beta.sdr=beta.sdr, d=d)
  return(ans)
}


######################
#for multivariate sir
######################
comp.sdr.prsirM = function( X, Y, d.true, H)
{ 
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  mn = 500
  
  M = matrix(0,p,p)
  Mhat = matrix(0,p,p)
  
  sigma = cov(X)
  eva = eigen(sigma)$values
  eve = eigen(sigma)$vectors
  sigmamrt = eve %*% diag(sqrt(1/eva)) %*% t(eve)
  sigmat = eve %*% diag(sqrt(eva)) %*% t(eve)
  
  xb = apply(X, 2, mean)
  xb = t(matrix(xb, p, n))
  
  x1 = X - xb
  z = as.matrix(x1)%*%sigmamrt        
  
  tpr = mvrnorm(mn, rep(0, q), diag(1, q))
  
  for (iter in 1:mn) 
  {
    t1 = tpr[iter, ]
    lt1 = sqrt((t(t1)%*%t1))
    t2 = t1/lt1
    ynew = Y%*%t2
    datanew = matrix(c(z,ynew),n,p+1) 
    datanew = datanew[order(datanew[,p+1]), ]
    
    znew = datanew[,1:p]
    nx = n/H
    zh = matrix(0, p, H)
    
    for (sl in 1:H)  
    {
      znewh = znew[((sl-1)*nx+1):(sl*nx), ]
      zh[ , sl] = colMeans(znewh)
    }
    
    Mhat = Mhat+zh%*%t(zh)/H
  }
  M = Mhat/mn
  
  chi.calc = rep(0,H)
  chi.tab = rep(0,H)
  for(i in 0:(H-2)){
    chi.calc[i+1]=n*sum(eigen(M)$val[(i+1):p])
    df=(H-i-1)*(p-i)
    chi.tab[i+1]=qchisq(0.05,df, ncp=0, lower.tail =  F, log.p = FALSE)}
  
  d=0
  k=1
  while(chi.calc[k]>chi.tab[k] & d < H-1){ d=d+1;k=k+1}
  
  if(d>d.true) d = d.true
  
  m = mat.sqrt(M)
  v = eigen(M)$vectors[,1:d]
  if(d == 1) v = matrix(v, ncol=1)
  beta.pprm<-sigmamrt %*% v
  beta.pprm<-apply(beta.pprm, 2, norm)
  ans = list(m=m, G=sigma, beta.sdr=beta.pprm, d=d)
  return(ans)
}


comp.sdr.euler<-function(X, y, d.true,n.w=1000,percent=0.2,s=0.02){
  M = alg(X,Y=y,s=0.02,n.w=1000,percent=0.2,method="euler")
  m = mat.sqrt(M)
  v = eigen(M)$vectors[,1:d.true]
  if(d.true == 1) v = matrix(v, ncol=1)
  beta.pprm<-v
  beta.pprm<-apply(beta.pprm, 2, norm)
  ans = list(m=m, G=var(X), beta.sdr=beta.pprm, d=d.true)
  return(ans)
}


comp.sdr.preuler<-function(X, y, d.true,n.w=1000,percent=0.2,s=0.02){   # X, y, d.true,n.w=1000,percent=0.2 
  M = alg(X,Y=y,s=0.02,n.w=1000,percent=0.2,method="euler")
  d <- permtest(X,Y=y,n.w=1000,percent=0.2,cutoff=0.05,way=c("euler"), alpha=0.05,ptime=100)
  m = mat.sqrt(M)
  v = eigen(M)$vectors[,1:d]
  if(d == 1) v = matrix(v, ncol=1)
  beta.pprm<-v
  beta.pprm<-apply(as.matrix(beta.pprm), 2, norm)
  ans = list(m=m, G=var(X), beta.sdr=beta.pprm, d=d)
  return(ans)
  
}

comp.sdr.kernel<-function(X, y, d.true){
  M = alg(X,Y=y,s=0.02,method="kernel")
  m = mat.sqrt(M)
  v = eigen(M)$vectors[,1:d.true]
  if(d.true == 1) v = matrix(v, ncol=1)
  beta.pprm<-v
  beta.pprm<-apply(beta.pprm, 2, norm)
  ans = list(m=m, G=var(X), beta.sdr=beta.pprm, d=d.true)
  return(ans)
}

comp.sdr.prkernel<-function(X, y, d.true){
  M = alg(X,Y=y,s=0.02,method="kernel")
  d <- permtest(X,Y=y,n.w=500,percent=0.2,cutoff=0.05,way=c("kernel"), alpha=0.05,ptime=100)  
  m = mat.sqrt(M)
  v = eigen(M)$vectors[,1:d]
  if(d == 1) v = matrix(v, ncol=1)
  beta.pprm<- v
  beta.pprm<-apply(as.matrix(beta.pprm), 2, norm)
  ans = list(m=m, G=var(X), beta.sdr=beta.pprm, d=d)
  return(ans)
}

######################################################################################
# functions for sparse sdr wrt G inner product							 #
######################################################################################

comp.sdr<-function(X, y, Xcat, d.true, method, ncat,nslices)  #comp.sdr(X, y, Xcat, d.true , method, ncat,nslices)
{
  out.m<-switch(method,
                sir    = comp.sdr.sir(X,y, nslices, d.true),
                prsirM = comp.sdr.prsirM(X,y, d.true, nslices),
                euler = comp.sdr.euler(X,y,d.true,n.w=1000,percent=0.2,s=0.02),
                preulerM = comp.sdr.preuler( X, y, d.true,n.w=1000,percent=0.2),
                kernel = comp.sdr.kernel( X, y, d.true),
                prkernelM = comp.sdr.prkernel( X, y, d.true)
  )
  ans<-list(m=out.m$m, G=out.m$G, beta.sdr=out.m$beta.sdr, d = out.m$d)
  return(ans)
}

# this function serves as a wrapper to call Hui Zou's solvebeta function
solve.beta<-function(x, y, G2, lambda1, lambda2)
{  
  # transform data to an L1 problem
  x.star<-rbind(x, sqrt(lambda2) * G2) 
  y.star<-c(y, rep(0, ncol(x)))
  # call solvebeta() with lambda2=0
  beta.est<-solvebeta(x.star, y.star, paras=c(0, lambda1), sparse="penalty")
  return(beta.est)
}


ssdr.lambda<-function(X, y, Xcat, ncat, method=c("sir","prsirM","euler","preulerM","kernel","prkernelM"), nslices=5, d.true,
                      lambda1=NULL, lambda2=1e-6, max.iter=200, eps.conv=1e-3)                
{
  n<-nrow(X)
  p<-ncol(X)
  method<-match.arg(method)
  # compute sdr components
  out.m<-comp.sdr(X, y, Xcat, d.true , method, ncat,nslices)
  d<-out.m$d
  m<-out.m$m
  M<-t(m) %*% m   
  G<-out.m$G
  G2<-mat.sqrt(G)
  G2.inv<-mat.sqrt.inv(G2)
  if(d!=0){
    # initial estimate of alpha and beta
    alpha<-out.m$beta.sdr      
    beta<-alpha
    for(i in 1:d) {
      ym<-m %*% alpha[,i]
      beta[,i]<-solve.beta(m, ym, G2, lambda1[i], lambda2)
    }
    # iteration
    iter<-0
    beta.n<-apply(beta, 2, norm)
    diff.conv<-1
    while((iter < max.iter) & (diff.conv[iter+1] > eps.conv)){
      z<-svd(G2.inv %*% M %*% beta)
      alpha<-G2.inv %*% (z$u) %*% t(z$v)
      for(i in 1:d) {
        ym<-m %*% alpha[,i]
        beta[,i]<-solve.beta(m, ym, G, lambda1[i], lambda2)
      }
      beta.n.new<-apply(beta, 2, norm)
      diff.conv<-c(diff.conv, max(abs(beta.n.new - beta.n)))
      beta.n<-beta.n.new
      iter<-iter + 1
    }
    # compute objective value
    comp1<-m %*% solve(G) - m %*% beta %*% t(beta)
    rss1<-sum(diag(comp1 %*% G %*% t(comp1)))
    
    z<-svd(G2.inv %*% M %*% beta)
    alpha<-G2.inv %*% (z$u) %*% t(z$v)
    comp2<-m %*% solve(G) - m %*% beta %*% t(alpha)
    rss2<-sum(diag(comp2 %*% G %*% t(comp2)))
    p.e<-sum(as.vector(beta) != 0)
    # normalize beta
    beta<-apply(beta, 2, norm)
    # return
    ans<-list(beta=beta, beta0=out.m$beta.sdr, alpha=alpha, diff.conv=diff.conv, iter=iter, p.e=p.e, rss1=rss1, rss2=rss2)
    ans1 = list(beta=beta,d=d)
    return(ans1)
  }
  else return(ans = list(d=d))
}


ssdr.wrap<-function(X, y, Xcat, ncat, method=c("sir","prsirM","euler","preulerM","kernel","prkernelM"), 
                    nslices=5, s1.range=NULL, s2.range=NULL, max.iter=200, eps.conv=1e-3)  
{
  n<-nrow(X)
  # compute criteria for all s1 and s2
  crit.all<-list()
  for(j in 1:length(s2.range)) {
    s2<-s2.range[j]
    
    crit<-NULL
    for(k in 1:length(s1.range)) {
      s1<-s1.range[k]
      
      out1<-ssdr.lambda(X, y, Xcat, ncat, method=method,nslices=nslices, lambda1=rep(s1, d), lambda2=s2)
      
      aic1<-n*out1$rss1 + 2 * out1$p.e
      bic1<-n*out1$rss1 + log(n) * out1$p.e
      out.crit<-c(aic1, bic1)
      crit<-cbind(crit, out.crit) 
    }
    crit.all[[j]]<-crit
  }
  
  # locate optimal s 
  s1.min<-NULL
  ct.min<-NULL
  for(j in 1:length(s2.range)) {
    s1.min.j<-ct.min.j<-NULL
    for(l in 1:length(out.crit)) {
      s1.min.j<-c(s1.min.j, s1.range[order(crit.all[[j]][l,])[1]])
      ct.min.j<-c(ct.min.j, min(crit.all[[j]][l,], na.rm=T))
    }
    s1.min<-rbind(s1.min, s1.min.j)
    ct.min<-rbind(ct.min, ct.min.j)
  }
  rownames(s1.min)<-as.character(s2.range); colnames(s1.min)<-c("aic1", "bic1")
  rownames(ct.min)<-as.character(s2.range); colnames(ct.min)<-c("aic1", "bic1")
  
  # beta estimate with given optimal s
  beta.est.all<-NULL
  s12.est.all<-NULL
  for(l in 1:length(out.crit)) {
    pos<-order(ct.min[, l])[1]
    s1<-s1.min[pos, l]
    s2<-s2.range[pos]
    out1<-ssdr.lambda(X, y, Xcat, ncat, method=method, d=d, nslices=nslices, lambda1=rep(s1, d), lambda2=s2)
    beta.est.all<-cbind(beta.est.all, out1$beta)
    s12.est.all <-cbind(s12.est.all,  c(s1, s2))
  }
  rownames(s12.est.all)<-c("s1", "s2"); colnames(s12.est.all)<-c("aic1", "bic1")
  
  ans<-list(beta.est.all=beta.est.all, beta.est0=out1$beta0, s12.est.all=s12.est.all, crit.all=crit.all, s1.min=s1.min, ct.min=ct.min)
  return(beta.est.all[,1:d])
}


########################################################################################
#	input : 												   #
#			X : design matrix (n x p)							   #
#			Y : the response vector (n x 1)						   #
#			d : the dimension of the cetral subspace					   #
#			H : the number of slices 							   #
#			m : the number of predictors considred at each step (p1 in the paper)#
####################################################
seq_sdr = function(X,Y,d,H,m, way=c("sir", "euler", "kernel")){    #X=newx;Y=y;d;H;m;way
  beta = diag(1,ncol(X))
  while	(	ncol(X)	>	m	)	
  {
    beta.step = c()
    step = floor(ncol(X)/m)
    newX = c()
    for(i in 1:step)   #i=1
    {
      step.i = min(i*m,ncol(X))
      X1 = X[,(1+(i-1)*m):step.i]
      if(step.i < ncol(X)) X2 = cbind(X[,(step.i+1):ncol(X)],newX)
      if(step.i == ncol(X)) X2 = newX
      
      if(i == step) X1 = X[,(1+(i-1)*m):ncol(X)]
      y = cbind(Y, X2)
      
      if (way=="sir"){
        prsirM.dr = ssdr.lambda(X1,y,Y, ncat=0, method='prsirM', nslices=H, d.true=d,
                                lambda1=c(1e-8,1e-8,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6),
                                lambda2=1e-6, max.iter=200, eps.conv=1e-3)
        } else if(way=="euler"){
        prsirM.dr = ssdr.lambda(X1,y,Y, ncat=0, method='preulerM', nslices=H, d.true=d,
                                lambda1=c(1e-8,1e-8,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6), 
                                lambda2=1e-6, max.iter=200, eps.conv=1e-3)
        } else {prsirM.dr = ssdr.lambda(X1,y,Y, ncat=0, method='prkernelM', nslices=H, d.true=d,
                                    lambda1=c(1e-8,1e-8,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6),
                                    lambda2=1e-6, max.iter=200, eps.conv=1e-3)}
      
      if(prsirM.dr$d==0)
      {
        beta.step.i = c()
        beta1.X = c()
      }
      else
      {
        beta1 = prsirM.dr$beta
        
        beta1.X = X1%*%beta1
        beta.step.i = matrix(0,ncol(X),ncol(beta1))
        if(i != step) beta.step.i[(1+(i-1)*m):step.i,] = beta1
        if(i == step) beta.step.i[(1+(i-1)*m):ncol(X),] = beta1
      }
      newX = cbind(newX,beta1.X)
      beta.step = cbind(beta.step,beta.step.i)
    }
    beta = beta%*%beta.step
    X = newX
  }
  if(ncol(newX)>d)
  {
    if(way=="sir"){
      prsirM.dr = ssdr.lambda(newX,y=Y,Y, ncat=ncat, method='sir', nslices=H, d.true = d,
                              lambda1=c(1e-8,1e-8,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6), 
                              lambda2=1e-6, max.iter=200, eps.conv=1e-3)}
    else if(way=="euler"){
      prsirM.dr = ssdr.lambda(newX,y=Y,Y, ncat=ncat, method='euler', nslices=H, d.true = d,
                              lambda1=c(1e-8,1e-8,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6),
                              lambda2=1e-6, max.iter=200, eps.conv=1e-3)}
    else{prsirM.dr = ssdr.lambda(newX,y=Y,Y, ncat=ncat, method='kernel', nslices=H, d.true = d,
                                 lambda1=c(1e-8,1e-8,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6),
                                 lambda2=1e-6, max.iter=200, eps.conv=1e-3)}
    beta = beta%*%prsirM.dr$beta
  }
  return(beta)
}


######################################################################################
# distance correlation between a random variable x and y
# function dist_corr accepts a vector of predictors X and 
# a scalar y and returns a vector with the distance correlations between the response 
# and the predictors
######################################################################################

dist_corr = function(X,y)
{
  n = nrow(X)
  dcor_xy = rep(0,ncol(X))
  
  for(k in 1:ncol(X))
  {
    x = X[,k]
    Xnew = Ynew = matrix(0,n,n)
    X_dup = matrix(x,n,n)
    Y_dup = matrix(y,n,n)
    for(i in 1:n){
      Xnew[,i] = x[i]
      Ynew[,i] = y[i]}
    
    X_dist = abs(X_dup - Xnew)
    Y_dist = abs(Y_dup - Ynew)
    a_xy = sum(X_dist*Y_dist)/n^2
    a_xx = sum(X_dist*X_dist)/n^2
    a_yy = sum(Y_dist*Y_dist)/n^2
    b_xy = sum(X_dist)/n^2 * sum(Y_dist)/n^2
    b_xx = sum(X_dist)/n^2 * sum(X_dist)/n^2
    b_yy = sum(Y_dist)/n^2 * sum(Y_dist)/n^2
    
    c_xy = c_xx = c_yy = 0
    for(i in 1:n){
      c_xy = c_xy + sum(X_dist[i,]*Y_dist)
      c_xx = c_xx + sum(X_dist[i,]*X_dist)
      c_yy = c_yy + sum(Y_dist[i,]*Y_dist)
    }
    c_xy = c_xy/n^3
    c_xx = c_xx/n^3
    c_yy = c_yy/n^3
    
    distcov_xy = a_xy + b_xy - 2*c_xy
    distcov_xx = a_xx + b_xx - 2*c_xx
    distcov_yy = a_yy + b_yy - 2*c_yy
    
    dcor_xy[k] = distcov_xy/sqrt((distcov_xx*distcov_yy))
  }
  return(dcor_xy)
}

########################################################################################
#	input : 												   #
#			x : design matrix (n x p)							   #
#			y : the response vector (n x 1)						   #
#			h : the number of slices used to order the predictors			   #
#			d : the dimension of the cetral subspace					   #
#			H : the number of slices used in the estimation 			   #
#			m : the number of predictors considred at each step (p1 in the paper)#
#	output : the basis for the central subspace (beta)					   #
########################################################################################

seq_sdr_with_ordering = function(x,y,h,d,H,m,way)    #h=4;d=1;H=10;m=20;x=X;y=y;way="kernel"                                
{
  con_cov = dist_corr(x,y)
  index = order(con_cov,decreasing=F)
  newx = x[,index[1:ncol(x)]]
  betaa = seq_sdr(newx,y,d,H,m,way)                                                    
  betahat = matrix(0,ncol(x),d)
  for(i in 1:ncol(x)){betahat[index[i],]= betaa[i,]}
  return(betahat)	
}

