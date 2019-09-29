# TODO: Add comment
# 
# Author: jh
###############################################################################

#' Internal fitting function
#'
#' @param y a multidimensional input vector.
#' @param g the number of groups.
#' @param error the desired relative error. (default 1e-5)
#' @param itermax the maximum of iterations. (default 1000)
#' @param pro1 the number of groups.
#' @param mu_neu the desired relative error. (default 1e-5)
#' @param Sigma_neu the maximum of iterations. (default 1000)
#' @param delta_neu of finding initial values. (default "kmeans")
#' @param nu_neu of finding initial values. (default "kmeans")
#' @param verbose of finding initial values. (default "kmeans")
#' @param mcmc calculate moments of truncated t distriubtion using mcmc algoirthm.
#' @param ... other inputs
#' @return a object of the class fitmixst.
#' @keywords fit mixed skew t.
#' @export

mmixstEst <- function(y, g, itermax, error, pro1, mu_neu, Sigma_neu, delta_neu, nu_neu, verbose=F, mcmc=F,...)
{
  p <- length(y[1,])
  n <- length(y[,1])
  lik_ges <- lik_ges_neu <- 1  
  
  #function for estimateing nu
  nuf <- function(nu,pro2,e2,e1) log(nu/2)-digamma(nu/2)-1/sum(pro2)*sum((pro2*(e2-e1-1))) 
  
  iter <- 0
  
  while(iter < itermax && abs(lik_ges / lik_ges_neu-1) > error || iter < 2)
  {
    mu <- mu_neu
    delta_neu1<-delta <- delta_neu
    Sigma <- Sigma_neu
    nu <- nu_neu
    Delta <- delta
    
    Omega <- list()
    Omega_inv <- list()
    Lambda <- list()
    q <- d <- e1 <- e2 <- t1 <- t2 <- t3 <- ew1 <- ew2 <- ew3 <- ew4 <- tmp <- S2 <-S3 <-pro2 <- matrix(NA,n,g)
    q <- y_star <- e3 <- ew3 <- array(NA,c(n,p,g))
    e4 <- ew4 <- array(NA,c(n,p,p,g))
    
    sumpro2 <- dmixst(y,list(pro=pro1,mu=mu,Sigma=Sigma,delta=delta,nu=nu))
    
    for(i in 1:g){
      
      Delta[[i]] <- diag(delta[[i]])
      Omega[[i]] <- Sigma[[i]]+(Delta[[i]])%*%t((Delta[[i]]))
      Omega_inv[[i]] <- solve(Omega[[i]])
      Lambda[[i]]<-diag(p)-t((Delta[[i]]))%*%Omega_inv[[i]]%*%(Delta[[i]])
      

      
      pro2[,i] <- dmixst(y,list(pro=pro1[i],mu=mu[i],Sigma=Sigma[i],delta=delta[i],nu=nu[i])) /sumpro2
      
      lik_ges <- lik_ges_neu
      lik_ges_neu <- sum(log(sumpro2))
      
      a<-array(NA,c(50,p,n))
      b<-array(NA,c(50,n))
      
      
      for(j in 1:n){
        
        q[j,,i]<-(Delta[[i]]) %*% (Omega_inv[[i]])%*%(y[j,]-mu[[i]])
        d[j,i]<-t(y[j,]-mu[[i]]) %*% (Omega_inv[[i]])%*%(y[j,]-mu[[i]])
        y_star[j,,i]<-q[j,,i] * sqrt((nu[[i]]+p) / (nu[[i]]+d[j,i]))
        
        
        if(!mcmc){
          moments <- truncatedt(a=rep(0,p), mu= -q[j,,i] , sigma= ((nu[[i]]+d[j,i])/(nu[[i]]+p + 2)) * Lambda[[i]], nu=round(nu[[i]]) + p +2)
          
          #e2[j,i]<- (nu[[i]] + p)/(nu[[i]] + d[j,i]) * 
          #  pmt(x=q[j,,i] * sqrt((nu[[i]]+p+2)/(nu[[i]]+d[j,i])), mean=rep(0,p), S=Lambda[[i]], df=round(nu[[i]]) + p + 2) / 
          #  pmt(x=y_star[j,,i], mean=rep(0,p), S=Lambda[[i]], df=round(nu[[i]]) + p)
          
          a <- (nu[[i]] + p)/(nu[[i]] + d[j,i])
          b <- pmt(x=q[j,,i] * sqrt((nu[[i]]+p+2)/(nu[[i]]+d[j,i])), mean=rep(0,p), S=Lambda[[i]], df=round(nu[[i]]) + p + 2)
          c <- pmt(x=y_star[j,,i], mean=rep(0,p), S=Lambda[[i]], df=round(nu[[i]]) + p)
          e2[j,i] <- a * b / c
          #cat(c, "\n")
          
          e1[j,i] <- e2[j,i] - log((nu[[i]] + d[j,i]) / 2) - (nu[[i]] + p) / (nu[[i]] + d[j,i]) +
            digamma((nu[[i]] + p) / 2)
          
          e3[j,,i] <-  - e2[j,i] * moments[[1]]
          e4[j,,,i] <- e2[j,i] * moments[[2]]
          
        }else{
          
          a[,,j] <- rtmvt(50,mean=q[j,,i],sigma=c(d[j,i]+nu[[i]])/(p+nu[[i]])*Lambda[[i]],
                          df=round(nu[[i]]+p),lower=rep(0,p),algorithm="gibbs")
          
          
          for(k in 1:50)
            b[k,j] <- (rgamma(1,shape=(nu[[i]]+2*p)/2,rate=c(t(a[k,,j]-q[j,,i])%*%solve(Lambda[[i]])%*%(a[k,,j]-q[j,,i])+d[j,i]+nu[[i]])/2))
          
          
          e2[j,i] <- mean(b[,j])
          e1[j,i] <- mean(log(b[,j]))
          
          bla1 <- matrix(NA,50,p)
          bla2 <- array(NA,c(50,p,p))
          
          for(k in 1:50)
          { 
            bla1[k,] <- c(a[k,,j] * b[k,j])  
            bla2[k,,] <- b[k,j] * a[k,,j] %*% t(a[k,,j])          
          }
          
          e3[j,,i] <- apply(bla1, 2, mean)
          e4[j,,,i] <- apply(bla2, 2:3, mean)
        }
        
        
      }
      
      
      #updating parameters
      
      mutmp = 0
      mutmp1 = 0
      for(k in 1:n){
        mutmp = mutmp + pro2[k,i]*e2[k,i] * y[k,] 
        mutmp1 = mutmp1 + pro2[k,i]*e3[k,,i]
      }
      
      mu_neu[[i]] = (mutmp - Delta[[i]] %*% mutmp1) / sum(pro2[,i]*e2[,i])

      
      deltmp1 = 0
      deltmp2 = 0
      for(k in 1:n){
        deltmp1 = deltmp1 +  pro2[k,i] * e3[k,,i] %*% t(y[k,] - mu_neu[[i]])
        deltmp2 = deltmp2 + pro2[k,i] * e4[k,,,i]
      }
      
      
      #new
      sigInv <- solve(Sigma_neu[[i]])
      
      
      #aus mclachlan 2012 (stimmt mit lin überein)
      #delta_neu[[i]] <- c(solve(sigInv * deltmp2) %*% diag(sigInv %*% deltmp1 ))
      
      #aus Lin 2010
      delta_neu[[i]] <- c(solve(sigInv * deltmp2) %*% (sigInv *deltmp1) %*% matrix(1, nrow=p, ncol=1))
      
      
      diagDel <- diag(delta_neu[[i]]) 
      
      #beide ok (von McLachlan 2012)
      #       sigtmp = 0
      #       for(k in 1:n){
      #         cenmu <- y[k,] - mu_neu[[i]]
      #         sigtmp  = sigtmp +  pro2[k,i] * ( diagDel %*% (e4[k,,,i]) %*% diagDel  - 
      #                                   (cenmu) %*% t(e3[k,,i]) %*% diagDel  -
      #                                   diagDel  %*%  (e3[k,,i]) %*% t(cenmu) +
      #                                   (cenmu) %*% t(cenmu)  * e2[k,i] )                                   
      #       }
      #       
      #        cat("sig ",sigtmp, "\n")
      
      #von Lin 2010
      sigtmp = 0
      for(k in 1:n){
        cenmu <- y[k,] - mu_neu[[i]]
        sigtmp  = sigtmp +  pro2[k,i] *(cenmu) %*% t(cenmu)  * e2[k,i] 
      }
      
      sigtmp = sigtmp + diagDel %*% deltmp2  %*%diagDel -
        diagDel %*% deltmp1 - t(deltmp1) %*% diagDel
      
      ind <- lower.tri(sigtmp)
      sigtmp[ind] <- t(sigtmp)[ind]
      
      
      Sigma_neu[[i]] <- 1 / sum(pro2[,i]) * sigtmp
      
      nu_neu[[i]]<-uniroot(function(nu) nuf(nu, pro2[,i], e2[,i], e1[,i]),interval=c(0.1,10e7))$root    
      
    }
    
    
    pro1<-apply(pro2,2,sum)/n
    
    iter <- iter+1
    
    #lik_ges <- lik_ges_neu
    
    #lik_ges_neu <- sum(log(dmixst(y,list(pro=pro1,mu=mu_neu,Sigma=Sigma_neu,delta=delta_neu,nu=nu_neu))))
    
    
    if(verbose) cat("Iteration: ", iter, "\nrelative error", abs(lik_ges / lik_ges_neu-1), "\nLoglikelihood: ", lik_ges_neu, 
                    "\npro", unlist(pro1), "\nmu", unlist(mu_neu), "\nSigma", unlist(Sigma_neu), "\ndelta", unlist(delta_neu), "\nnu", unlist(nu_neu), "\n\n")
  }
  
  
  
  
  #sort by size for easy comparison of models
  # 	l <- order(mu_neu, decreasing=T)
  # 	mu_neu <- mu_neu[l]
  # 	delta_neu <- delta_neu[l]
  # 	Sigma_neu <- Sigma_neu[l]
  # 	pro1 <- pro1[l]
  # 	pro2 <- pro2[,l]
  
  list(mu=mu_neu, Sigma=Sigma_neu, delta=delta_neu, nu=nu_neu, pro=pro1, iter=iter, logLik=lik_ges_neu, 
       posteriori=pro2, empcov = matrix(NA,0,0), g=g, p=p)
}

