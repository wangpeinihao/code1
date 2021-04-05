library(MASS)
source("algorithm.R")
source("criteria.R")


######################
#generate model
#####################
n <- 200
p <- 10
mu <- rep(0,p)
sigma <- diag(1,p)
x <- mvrnorm(n, mu, sigma)
beta <- c(1, 1, rep(0,p-2))
epsilon <- rnorm(n)
y <- (x%*%beta)^2 + 0.4*epsilon                                    
sbeta <- svd(beta)$u  # standardize the original direction
d <- 1

############################
# SDR example
###########################
#euler
M <- alg(X=x,Y=y,n.w=20000,percent=0.2,method="euler")  
eig <- eigen(M)   
eval <- eig$values       
evec <- eig$vectors      
odir <- evec[,c(1:d)]   
dir <- svd(odir)$u       # standardize the estimated direction.
euler.res <- crifunc(d=d,dir=dir,sbeta = sbeta)

#kernel
M <- alg(X=x,Y=y,method="kernel")
eig <- eigen(M)   
eval <- eig$values       
evec <- eig$vectors      
odir <- evec[,c(1:d)]   
dir <- svd(odir)$u      
kernel.res <- crifunc(d=d,dir=dir,sbeta = sbeta)


###################################
# Alasso example
###################################
source("alasso.R")
kernelest =  al.ssdr(method="kernel",x=x,y=y, d, H, lambda1.range= seq(0.0001,0.01,0.0001),lambda2.range,delta=0.001)[[2]] 
eulerest =  al.ssdr(method="euler",x=x,y=y, d, H,  lambda1.range = seq(0.001,1,0.03),lambda2.range,delta=0.001)[[2]]

#############################
#CISE example
#############################
source("cise.R")
kernelest0 =  cise(y,x,di=d,pw= seq(0.0001,0.5,0.001),method="KERNEL")
kernelest =  kernelest0$st
eulerest0 =  cise(y,x,di=d,pw= seq(0.001,2,0.03),method="EULER") 
eulerest=eulerest0$st

###########################
# permutation test example
###########################
source("permtest.R")
estdkernel <- permtest(X=x,Y=y,n.w=1000,percent=0.2,cutoff=0.05,way="kernel", alpha=0.05,ptime=100)
estdeuler <- permtest(X=x,Y=y,n.w=1000,percent=0.2,cutoff=0.05,way="euler", alpha=0.05,ptime=100)

########################################################
# variable screening example
# R function "DC.R" is requsted from Yang et al.(2019)
# Use the two-stage selection procedure in DC.R
#######################################################
source("DC.R")  
n = 200;p = 1000;d = 1 
active = c(100,200,300,400,500,600,700,800)
inactive = c(seq(1,99),seq(101,199),seq(201,299),seq(301,399),seq(401,499),seq(501,599),seq(601,699), seq(701,799),seq(801,1000))
beta = rep(0,p); beta[active]= sqrt(7/8); beta[inactive]= 0
sbeta <- beta/(sqrt(sum(beta^2)))
Y<- matrix(rnorm(n,0,0.5),n,1)
epsilon <- mvrnorm(n, rep(0,p), diag(1,p))           
X1 <- Y%*%t(beta)+0.5*epsilon

TH<-VarSelect_two(Y,X1,floor(n/log(n)))$DC[,1]    #DC-Two
locrep<-TH[order(TH)]
X <- X1[,locrep]  #dim(X)
beta.x = X1%*%sbeta

out.euler <- alg(X=X,Y=Y,n.w=20000,percent=0.2,method="euler")  
betahate <- svd(out.euler)$u[,1:d]
edir <- matrix(0,p,d)
edir[locrep,] <- betahate 
correx = abs(cor(beta.x,betae.x))
correbt <- crifunc(d, dir=as.matrix(edir), sbeta)
euler.res <- c(correx,correbt)


out.kernel <- alg(X=X,Y=Y,n.w=1,percent=1,method="kernel")
betahatk <- svd(out.kernel)$u[,1:d]
betak.x <- X%*%betahatk     
kdir <- matrix(0,p,d)
kdir[locrep,] <- betahatk 
corrkx = abs(cor(beta.x,betak.x))
corrkbt <- crifunc(d, dir=as.matrix(kdir), sbeta)
kernel.res <- c(corrkx,corrkbt)


###################################################
# seqential SVS example
###################################################
source("seqentialSVS.R")
library(Rlab)
library(dr)
library(lars)
n = 200;p = 1000
active = c(100,200,300,400,500,600,700,800)
inactive = c(seq(1,99),seq(101,199),seq(201,299),seq(301,399),seq(401,499),seq(501,599),seq(601,699),
             seq(701,799),seq(801,1000))
beta = rep(0,p);beta[active]= sqrt(7/8);beta[inactive]= 0
tpr = rep(0,length(active))
fpr = rep(0,length(inactive))
mu = rep(0,p)
sigma = diag(1,p)
X = mvrnorm(n, mu, sigma)
y = X%*%beta + 1.5*rnorm(n,0,1)
betahat =  seq_sdr_with_ordering(X,y,4,1,10,20,way="euler")
#betahat =  seq_sdr_with_ordering(X,y,4,1,10,20,way="kernel")
beta.x = X%*%beta
betahat.x =X%*%betahat
corr = cor(beta.x,betahat.x)
A = betahat%*%solve(t(betahat)%*%betahat)%*%t(betahat)-beta%*%solve(t(beta)%*%beta)%*%t(beta)
deltaf = sqrt(sum(diag(A%*%t(A))))
deltam = max(abs(eigen(A)$value))
for(k in 1:length(active)){if(betahat[active[k]]!=0) tpr[k]=1}
for(k in 1:length(inactive)){if(betahat[inactive[k]]!=0) fpr[k]=1}
result <- c(corr,deltaf,deltam,tpr,fpr)





