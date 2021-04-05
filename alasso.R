#######################################################
#function used to get the adaptive lasso sparse estimation
#we modified the codes provided by Wu and Yin(2015)
#For the sir, save, phd method, using Wu and Yin's code
#Wu and Yin's code is a modification of Li(2007)
#have to have packages glmnet/elasticnet downloaded
######################################################
library(glmnet)
library(elasticnet)
source("algorithm.R")

#####################################
# center a matrix
#####################################
center = function(x){
  return(t(t(x)-apply(x,2,mean)))
}
###################################
#standardize a matrix
###################################
standardize = function(x){
  return(center(x)%*%mppower(cov(x),-0.5,0.000000001))
}

standardize_col = function(A){
  for (i in 1:ncol(A)){
    if (sqrt(drop((t(A[,i])%*%A[,i])))!=0) A[,i]= A[,i]/sqrt(drop((t(A[,i])%*%A[,i])))
  }
  return(A)
}


switch_col = function(A, i,j){
  temp = A[,i]
  A[,i]=A[,j]
  A[,j]=temp
  return(A)
}

sign_beta = function(a,b){
  a <-as.matrix(a)
  b <-as.matrix(b)
  for (i in 1:ncol(b)){
    if (sum(abs(a[,i]-b[,i]))>sum(abs(-a[,i]-b[,i]))) b[,i]=-b[,i]
  }
  return(b)
}

order_beta = function(a,b){
  a <-as.matrix(a)
  b <-as.matrix(b)
  p = nrow(b)
  d = ncol(b)
  for (i in 1:d){
    v = numeric(d)
    for (j in i:d){
      v[j]=vc(matrix(a[,i],p,1),matrix(b[,j],p,1))
    }
    b = switch_col(b,i,min(which(v==max(v))))
  }
  return(b)
}


### compute vector correlation between matrix beta and beta_hat ###
vc = function(beta, beta_hat){
  beta = standardize_col(beta)
  beta_hat = standardize_col(beta_hat)
  return(sqrt(abs(det(t(beta)%*%beta_hat%*%t(beta_hat)%*%beta))))
}

### obtaining S.I.R. kernel matrix by evenly split the response range ###
get_sir_matrix_even_range = function(x,y,H){
  z = standardize(x)
  n = length(y)
  p = ncol(x)
  sliced_y = numeric(n)
  slice_prop = numeric(H)
  divpt = c(min(y)-0.1,min(y)+((range(y)[2]-range(y)[1])/H)*c(1:H))
  sir_matrix = matrix(0,p,p)
  for(i in 1:H){
    sliced_y[(y<=divpt[i+1])&(y>divpt[i])]=i
    slice_prop[i] = sum(sliced_y[sliced_y==i])/(n*i)
    slice_mean = matrix(colMeans(matrix(z[sliced_y==i,],length(z[sliced_y==i,])/p,p)))
    if (slice_prop[i]!=0)
      sir_matrix = sir_matrix + slice_prop[i]*slice_mean%*%t(slice_mean)
  }
  return(sir_matrix)
}
### obtaining S.I.R. estimate for a given dimension ###
sir_even_range = function(x,y,H,d){
  p = ncol(x)
  M = get_sir_matrix_even_range(x,y,H)
  beta_hat=mppower(cov(x),-0.5,0.000000001)%*%(eigen(M)$vectors)[,1:d]
  beta_hat = standardize_col(matrix(beta_hat,p,d))
  return(beta_hat)
}

### Principal Hessian Directions (P.H.D.) ####
### obtaining P.H.D. using residuals ###
get_phd_matrix_r = function(x,y){
  r = residuals(lm(y~x))
  xbar=colMeans(x)
  n = nrow(x)
  p = ncol(x)
  sigma_rxx = matrix(0,p,p)
  for(i in 1:n){
    sigma_rxx = sigma_rxx +r[i]*t(t(x[i,]-xbar))%*%t(x[i,]-xbar)
  }
  m_phd = solve(cov(x))%*%(sigma_rxx/n)
  
  return(m_phd)
}
get_phd_matrix = function(x,y){
  r = residuals(lm(y~x))
  xbar=colMeans(x)
  n = nrow(x)
  p = ncol(x)
  sigma_rxx = matrix(0,p,p)
  for(i in 1:n){
    sigma_rxx = sigma_rxx +r[i]*t(t(x[i,]-xbar))%*%t(x[i,]-xbar)
  }
  m_phd = solve(cov(x))%*%(sigma_rxx/n)
  
  return(m_phd%*%t(m_phd))
}
### obtaining P.H.D. estimate for a given dimension ###
phd_r = function(x,y,d){
  p = ncol(x)
  m_phd = get_phd_matrix_r(x,y)
  result = mppower(cov(x),-0.5,0.000000001)%*%eigen(m_phd)$vectors[,1:min(d,p)]
  return(standardize_col(matrix(result,p,d)))
}


### obtaining S.A.V.E. kernel matrix by evenly split the response range ###
get_save_matrix_even_range = function(x,y,H){
  z = standardize(x)
  n = length(y)
  p = ncol(x)
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
  return(save_matrix)
}
### obtaining S.A.V.E. estimate for a given dimension ###
save_even_range = function(x,y,H,d){
  p = ncol(x)
  M = get_save_matrix_even_range(x,y,H)
  beta_hat=mppower(cov(x),-0.5,0.000000001)%*%(eigen(M)$vectors)[,1:d]
  beta_hat = standardize_col(matrix(beta_hat,p,d))
  return(beta_hat)
}

kernel_est <- function(x,y,d,p){
  kernelres0 <- alg(X=x,Y=y,n.w=1000,method="kernel")
  eig <- eigen(kernelres0)   
  eval <- eig$values       
  evec <- eig$vectors
  odir <- evec 
  kernelres <- svd(as.matrix(odir)[,1:d])$u
  return(kernelres)
}


euler_est <- function(x,y,d,p){
  eulerres0 <- alg(X=x,Y=y,n.w=1000,percent=0.2,method="euler")
  eig <- eigen(eulerres0)   
  eval <- eig$values       
  evec <- eig$vectors      
  odir <- evec   
  eulerres <- svd(as.matrix(odir)[,1:d])$u
  return(eulerres)
}

kernel_M <- function(x,y){
  kernelM <- alg(X=x,Y=y,n.w=1000,method="kernel")
  return(kernelM)
}

euler_M <- function(x,y,s=0.02,n.w=1000,percent=0.2){
  eulerM <- alg(X=x,Y=y,n.w=1000,percent=0.2,method="euler")
  return(eulerM)
}

### sparse basis direction estimations ####
sparse = function(method,x, y, d, H, lambda1,lambda2,delta=0){ #out = sparse(method,x, y, d, H, w, lambda1,lambda2,delta)
  ### parameters ###
  n<-nrow(x)
  p<-ncol(x)
  G = cov(x)
  M = switch(	method,
              sir=get_sir_matrix_even_range(x,y,H), 
              phd=get_phd_matrix(x,y),
              save=get_save_matrix_even_range(x,y,H),
              kernel=kernel_M(x,y),
              euler=euler_M(x,y)
  ) + delta*G
  m = mppower(M,0.5,0.000000001) 
  alpha = switch(	method,
                  sir=sir_even_range(x,y,H,d), 
                  phd=phd_r(x,y,d),
                  save=save_even_range(x,y,H,d),     ##alpha value for iteration
                  kernel=kernel_est(x,y,d,p),
                  euler=euler_est(x,y,d,p)
  )
  alpha <- as.matrix(alpha)
  beta.est = alpha
  last = beta.est
  iter = 0
  diff = 1
  while (iter<100&diff>1e-3){
    for(i in 1:d) {
      u = m%*%alpha[,i]
      x.star=rbind(m, sqrt(lambda2) * mppower(G,-0.5,0.000000001)) 
      y.star=c(u, rep(0, ncol(x)))
      w = abs(lm(y.star~x.star)$coefficients[-1])
      x.star=scale(x.star,center=FALSE,scale=1/w)
      # call solvebeta() with lambda2=0
      beta.est[,i] = solvebeta(x.star, y.star, para=c(0, lambda1), sparse="penalty")
      beta.est[,i]=beta.est[,i]*w
    }
    u = svd(mppower(G,-0.5,0.000000001)%*%M%*%beta.est)$u
    v = svd(mppower(G,-0.5,0.000000001)%*%M%*%beta.est)$v
    alpha = mppower(G,-0.5,0.000000001)%*%u%*%t(v)
    iter = iter+1
    diff = max(abs(beta.est - last))
    last = beta.est
  }
  #beta.est=apply(beta.est, 2, norm)
  beta.est = standardize_col(beta.est)
  # compute objective value
  comp1<-m %*% solve(G) - m %*% beta.est %*% t(beta.est)
  rss1<-sum(diag(comp1 %*% G %*% t(comp1)))
  z<-svd(mppower(G,-0.5,0.000000001) %*% M %*% beta.est)
  alpha<-mppower(G,-0.5,0.000000001) %*% (z$u) %*% t(z$v)
  comp2<-m %*% solve(G) - m %*% beta.est %*% t(alpha)
  rss2<-sum(diag(comp2 %*% G %*% t(comp2)))
  p.e<-sum(as.vector(beta.est) != 0)
  return(list(beta.est=beta.est,rss1=rss1,rss2=rss2,p.e=p.e))
}

### selecting optimal tuning parameter for sparse basis direction estimations ####
choose_tuning = function(method,x, y, d, H,lambda1.range,lambda2.range,delta=0){
  crit = matrix(0,length(lambda1.range)*length(lambda2.range),4)
  n = nrow(x)
  rn = 0
  for(j in 1:length(lambda2.range)) {
    lambda2=lambda2.range[j]
    for(k in 1:length(lambda1.range)) {
      lambda1=lambda1.range[k] 
      rn = rn + 1
      out = sparse(method,x, y, d, H,lambda1,lambda2,delta)  #w=rep(1,p) x==subx
      aic = n*out$rss1 + 2 * out$p.e
      bic = n*out$rss1 + log(n) * out$p.e
      crit[rn,1] = lambda1
      crit[rn,2] = lambda2
      crit[rn,3] = aic
      crit[rn,4] = bic
    }
  }
  lambda_aic = crit[crit[,3]==min(crit[,3]),]
  lambda_bic = crit[crit[,4]==min(crit[,4]),]
  return(lambda_bic)
}

###sparse basis direction estimations ####
al.ssdr = function(method,x,y,d, H=6,lambda1.range,lambda2.range,delta){
  n = nrow(x)
  p = ncol(x)
  ini.beta0 = switch(	method,
                      sir=sir_even_range(x,y,H,d), 
                      phd=phd_r(x,y,d),
                      save=save_even_range(x,y,H,d),
                      kernel=kernel_est(x,y,d,p),
                      euler=euler_est(x,y,d,p)
  )
  opt0 = choose_tuning(method,x, y, d, H, lambda1.range,lambda2.range=0,delta=0.01)
  sparse.beta0 = sparse(method,x, y, d, H, opt0[1],opt0[2],delta)$beta.est
  sparse.beta0 = order_beta(ini.beta0,sparse.beta0)
  sparse.beta0 = sign_beta(ini.beta0,sparse.beta0)
  return(list(ini.beta=ini.beta0,sparse.beta=sparse.beta0))
}


