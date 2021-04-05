##########################################################
# the function crifun is used to calculate the criteria
##########################################################
crifunc <- function(d, dir, sbeta){
  if(d == 1) {
    r <- sqrt(t(dir)%*%sbeta%*%t(sbeta)%*%dir)   #trace correlation, larger, better
  } else if(d != 1 ) {
    r <- sqrt(mean(svd(t(dir)%*%sbeta%*%t(sbeta)%*%dir)$d))
  }
  deltam <- max(svd(dir%*%t(dir)-sbeta%*%t(sbeta))$d) # maximum singular value, smaller, better
  deltaf<-sqrt(sum(diag((dir%*%t(dir)-sbeta%*%t(sbeta))%*%t(dir%*%t(dir)-sbeta%*%t(sbeta)))))  #Frobenius norm, smaller, better
  mi <- c()
  for(i in 1:d){
    mi[i] <-  sqrt(t(dir[,i]-sbeta%*%t(sbeta)%*%dir[,i])%*%(dir[,i]-sbeta%*%t(sbeta)%*%dir[,i])) #smaller, better
  }
  crit <- c(r,deltam,deltaf,mi)
  return(crit)
}
