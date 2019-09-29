# TODO: Add comment
# 
# Author: jh
###############################################################################

#' Internal function for moments of truncated t
#'
#' @param a truncation at
#' @param mu mean of t distribution
#' @param sigma scale of t distribution
#' @param nu degrees of freedom
#' @return first and second moment of the truncated t distribution
#' @keywords truncated t distribution
#' @export
#' 
truncatedt <- function(a,mu,sigma,nu){
  p <- length(a)
  m1 <- vector(length=p)
  h  <- m3 <- matrix(0,p,p)
  tmp <- 1
  tmp2 <- vector(length=p)
  
  for (i in 1:p){
    #evtl falsch in ohagan
    astar  <-  (a[-i]-mu[-i])-(a[i]-mu[i])*1/sigma[i,i]*sigma[,i][-i];
    sigmastar <- (nu+1/sigma[i,i]*(a[i]-mu[i])^2)/(nu-1)*(sigma[-i,-i]-1/sigma[i,i]*
                                                            sigma[,i][-i]%*%t(sigma[,i][-i]))
    m1[i] <-(2*pi*sigma[i,i])^(-0.5) * (nu/(nu+1/sigma[i,i]*(a[i]-mu[i])^2))^((nu-1)/2)*
      gamma((nu-1)/2)/gamma(nu/2)*
      sqrt(nu/2)*pmt(astar,mean=rep(0,p-1),S=sigmastar,nu-1)/pmt((a-mu),mean=rep(0,p),S=sigma,df=nu)
    
    for(j in 1:p){
      if(i!=j){
        
        #sigma
        nustar  <- nu+c(t(a[c(i,j)]-mu[c(i,j)])%*%solve(sigma[c(i,j),c(i,j)])%*%(a[c(i,j)]-mu[c(i,j)]))
        astarstar <- (a[c(-i,-j)]-mu[c(-i,-j)]) - cbind(sigma[,i][c(-i,-j)],sigma[,j][c(-i,-j)])%*%
          solve(sigma[c(i,j),c(i,j)])%*%(a[c(i,j)]-mu[c(i,j)])        
        if(p>2){
          sigmastarstar <- nustar/(nu-2)*(sigma[c(-i,-j),c(-i,-j)]-cbind(sigma[,i][c(-i,-j)],sigma[,j][c(-i,-j)])%*%
                                            solve(sigma[c(i,j),c(i,j)])%*%t(cbind(sigma[,i][c(-i,-j)],sigma[,j][c(-i,-j)])))
          tmp <- pmt(astarstar,mean=rep(0,p-2),sigmastarstar,nu-2)
        }
        
        
        h[i,j] <- -1/(2*pi*sqrt(sigma[i,i]*sigma[j,j]-sigma[i,j]^2))*(nu/(nu-2))*(nu/nustar)^(nu/2-1)*tmp/
          pmt(a-mu,mean=rep(0,p),S=sigma,df=nu)
        tmp2[i] <- tmp2[i] + sigma[i,j]*h[i,j]
      }  
    }
  }
  
  for(i in 1:p){
    h[i,i] <- 1/sigma[i,i]*(a[i]-mu[i])*m1[i]-1/sigma[i,i]*tmp2[i]
  }
  
  m2 <- mu-sigma%*%m1
  m3  <- -mu %*% t(mu) + mu %*% t(m2) + m2 %*% t(mu) + nu/(nu-2)*pmt(a-mu,mean=rep(0,p),S=nu/(nu-2)*
                                                                       sigma,df=nu-2)/pmt(a-mu,mean=rep(0,p),S=sigma,df=nu)*sigma -sigma%*%h%*%sigma
  list(moment1 = m2, moment2 = m3)
}
