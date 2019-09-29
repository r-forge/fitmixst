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

# library(fitmixst4)
# library(mnormt)
# y <- cbind(rnorm(10,0), rnorm(10,.5), rnorm(10,1))
# g <- 1
# d1 <- dim(y)
# itermax <- 1e3
# 
# #if(method =="kmeans"){
# 
# a <- 0.5
# init <- kmeans(y, g,nstart=5)
# 
# pro1 <- init$size/length(y[,1])
# mu_neu <- delta_neu <- Sigma_neu <- nu_neu <- list()
# #nu_neu <- vector()
# length(mu_neu) <- g
# length(delta_neu) <- g
# length(Sigma_neu) <- g
# length(nu_neu) <- g
# 
# for (j in 1:g) {
#   if(sum(init$cluster == j) > 1){
#     Sigma_neu[[j]] <- a*var(y[init$cluster == j,])
#     mu_neu[[j]] <- init$centers[j, ]
#     delta_neu[[j]] <- sign(sum((y[init$cluster == j,] -(rep(mu_neu[[j]], sum(init$cluster == j))))^3)) * 
#       sqrt((1-a)* diag(var(y[init$cluster == j,]))*(1-2/pi))
#     mu_neu[[j]] <- init$centers[j, ] - sqrt(2/pi) * delta_neu[[j]]
#     dimnames(Sigma_neu[j]) <- NULL
#     names(mu_neu[j]) <- NULL
#     names(delta_neu[j]) <- NULL
#     nu_neu[[j]] <- 40
#   }else{
#     Sigma_neu[[j]] <- diag(d1[2])
#     mu_neu[[j]] <- init$centers[j, ]
#     delta_neu[[j]] <- matrix(0, nrow=d1[2], ncol =1)
#     mu_neu[[j]] <- init$centers[j, ]
#     dimnames(Sigma_neu[j]) <- NULL
#     names(mu_neu[j]) <- NULL
#     names(delta_neu[j]) <- NULL
#     nu_neu[[j]] <- 40
#     warning("Parameter g to large. Too many groups specified!")
#   }
# }  
# 
# #}
# 
# i <- 1
#function for estimateing nu


mmixstEst2 <- function(y, g, itermax, error, pro1, mu_neu, Sigma_neu, delta_neu, nu_neu, verbose=F, mcmc=F,...)
{

  MlistVec <- function(Mlist, vec)
  {
    p <- nrow(Mlist[[1]])
    scal <- matrix(rep(vec, each=p*p), ncol=p, byrow=T)
    M2 <- matrix(unlist(Mlist), ncol=p, byrow=T)
    
    M2 * scal
  }
  
  blockSum <- function(longMat, n)
  {
    t(rep(1,n) %x% diag(ncol(longMat))) %*% longMat
  }
  
  # .vec und vec2 sind gleich in Sachen speed. daher vec2, da weniger
  # speicher gebraucht wird. Ab .vec3 machen sie schmarrn.
  pmt.vec <- function(mat, ...)
  {
    apply(mat, 1, pmt, ...)
  }
  pmt.vec2 <- function(mat, ...)
  {
    n <- nrow(mat)
    res <- double(n)
    for(i in 1:n) res[i] <- pmt(mat[i,],...)
    res
  }
  pmt.vec3 <- function(mat, ...)
  {
    lapply(t(mat), pmt, ...)
  }
  pmt.vec4 <- function(mat, ...)
  {
    sapply(t(mat), pmt, ...)
  }
  pmt.vec5 <- function(mat, ...)
  {
    vapply(t(mat), pmt, FUN.VALUE=1.22, ...)
  }
  
  
  p <- length(y[1,])
  n <- length(y[,1])
  lik_ges <- lik_ges_neu <- 1  
  
  
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
    
    for(i in 1:g)
      {
      
      Delta[[i]] <- diag(delta[[i]])
      Omega[[i]] <- Sigma[[i]]+(Delta[[i]])%*%t((Delta[[i]]))
      Omega_inv[[i]] <- solve(Omega[[i]])
      Lambda[[i]]<-diag(p)-t((Delta[[i]]))%*%Omega_inv[[i]]%*%(Delta[[i]])
      
      
      
      pro2[,i] <- dmixst(y,list(pro=pro1[i],mu=mu[i],Sigma=Sigma[i],delta=delta[i],nu=nu[i])) /sumpro2
      
      lik_ges <- lik_ges_neu
      lik_ges_neu <- sum(log(sumpro2))
      
      a<-array(NA,c(50,p,n))
      b<-array(NA,c(50,n))
      
      
      #for(j in 1:n){
      
      # orig
      #q[j,,i]<- (Delta[[i]]) %*% (Omega_inv[[i]])%*%(y[j,]-mu[[i]])
      # neu
      Ycenter <- t(t(y) - c(mu[[i]]))
      Yscaled <- Ycenter %*%  Omega_inv[[i]] 
      q[,,i] <-   Yscaled %*% Delta[[i]]
      #orig
      #d[j,i]<-t(y[j,]-mu[[i]]) %*% (Omega_inv[[i]])%*%(y[j,]-mu[[i]])
      # neu
      d[,i] <-  diag(Ycenter %*% t(Yscaled))
      # orig
      #y_star[j,,i]<-q[j,,i] * sqrt((nu[[i]]+p) / (nu[[i]]+d[j,i]))
      # neu
      y_star[,,i] <- q[,,i] * sqrt((nu[[i]]+p) / (nu[[i]]+d[,i]))
      
      
      
      if(!mcmc){
        #neu: hilfobjekte
        cat(dput(q[,,i]))
        ttmu <- lapply(seq_len(n), function(ii) -q[,,i][ii,])
        ttsigma <- lapply((nu[[i]]+d[,i])/(nu[[i]]+p + 2), "*", Lambda[[i]])
        
        # orig (#hängt nicht von j ab??)
        #          moments <- truncatedt(a=rep(0,p), 
        #                                 mu= -q[j,,i] , 
        #                                 sigma= ((nu[[i]]+d[j,i])/(nu[[i]]+p + 2)) * Lambda[[i]], 
        #                                  nu=round(nu[[i]]) + p +2)
        #neua als liste
        momentsList <- vector(length=n, mode="list")
        for(j in 1:n)
        {
          momentsList[[j]] <- truncatedt(
            sigma = ttsigma[[j]],
            mu = ttmu[[j]],
            a = rep(0,p),
            nu = round(nu[[i]]) + p +2)
        }
        
                           
        #momentsList <- mapply(truncatedt, 
        #                      sigma = ttsigma,
        #                      mu = ttmu,
        #                      MoreArgs = list(
        #                        a = rep(0,p),
        #                        nu = round(nu[[i]]) + p +2),
        #                      SIMPLIFY = FALSE
        #)
        #cat(sapply(momentsList, function(x) x[[2]]), "\n")
        #orig
       # for(j in 1:n)
       # {
      #    cat(q[j,,i] * sqrt((nu[[i]]+p+2)/(nu[[i]]+d[j,i])), "\n", "\n")
      #  }
        #e2[j,i]<-
        #           (nu[[i]] + p)/(nu[[i]] + d[j,i]) * 
        #             pmt(x=q[j,,i] * sqrt((nu[[i]]+p+2)/(nu[[i]]+d[j,i])), mean=rep(0,p), S=Lambda[[i]], df=round(nu[[i]]) + p + 2) / 
        #             pmt(x=y_star[j,,i], mean=rep(0,p), S=Lambda[[i]], df=round(nu[[i]]) + p)
        #neu
      #  pmtx <- lapply(seq_len(n), function(ii) 
       #   (diag(sqrt((nu[[i]]+p+2)/(nu[[i]]+d[,i]))) %*% q[,,i])[ii,])
        #cat(pmtx[[1]], "\n", "\n")
        bla <- rep(sqrt((nu[[i]]+p+2)/(nu[[i]]+d[,i])), each=p)
        bla2 <- matrix(bla, ncol=p, byrow=T)
        #cat(dim(bla2), "\n")
        #cat(dim(q[,,i]), "\n")
        pmtx <- (q[,,i]) * bla2
        #cat(pmtx2, "\n")
        #cat("Hallo", matrix(bla, ncol=p, byrow=T)[1,], "\n")
        #cat("Jettz", matrix(bla, ncol=p, byrow=T) * matrix(q[,,i], ncol=p, byrow=F), "\n")
        #cat(dim(pmtx), "\n")
        
        
        #pmtx2 <- lapply(seq_len(n), function(ii) 
        #  y_star[,,i][ii,])
        pmtx2 <- y_star[,,i]
        #cat(dim(pmtx2), "\n")
        
#         e2[,i] <- (nu[[i]] + p)/(nu[[i]] + d[,i]) * pro2[,i] *
#           (apply(pmtx, 2, pmt,
#                  mean=rep(0,p), 
#                  S=Lambda[[i]], 
#                  df=round(nu[[i]]) + p + 2)) /
#           (apply(pmtx2, 2, pmt, 
#                  mean=rep(0,p), S=Lambda[[i]], 
#                  df=round(nu[[i]]) + p))
        a <- (nu[[i]] + p)/(nu[[i]] + d[,i]) * pro2[,i]
        b <- pmt.vec2(mat = pmtx, 
                    mean=rep(0,p), 
                    S=Lambda[[i]], 
                    df=round(nu[[i]]) + p + 2)
        c <- pmt.vec2(mat = pmtx2, 
                    mean=rep(0,p), S=Lambda[[i]], 
                    df=round(nu[[i]]) + p)
        
        e2[,i] <- a* b / c
        #print(c)
        
        
        # orig
        #e1[j,i] <- e2[j,i] - log((nu[[i]] + d[j,i]) / 2) - (nu[[i]] + p) / (nu[[i]] + d[j,i]) +
        #  digamma((nu[[i]] + p) / 2)
        # neu
        e1[,i] <- (e2[,i] - log((nu[[i]] + d[,i]) / 2) - (nu[[i]] + p) / (nu[[i]] + d[,i]) +
          digamma((nu[[i]] + p) / 2) ) * pro2[,i]
        
        #orig  
        #  e3[j,,i] <-  - e2[j,i] * moments[[1]]
        #  e4[j,,,i] <- e2[j,i] * moments[[2]]
        #neu
        e3[,,i] <-  -t(sapply(momentsList, function(x) x[[1]])) * e2[,i]
        m2list <- lapply(momentsList, function(x) x[[2]])
        e4list <- MlistVec(Mlist=m2list,vec=e2[,i])
          #mapply("*", x=m2list, y=as.list(e2[,i]), SIMPLIFY = FALSE)
        # e4 nur als Liste !!
        
        
      }else{
        #           
        #           a[,,j] <- rtmvt(50,mean=q[j,,i],sigma=c(d[j,i]+nu[[i]])/(p+nu[[i]])*Lambda[[i]],
        #                           df=round(nu[[i]]+p),lower=rep(0,p),algorithm="gibbs")
        #           
        #           
        #           for(k in 1:50)
        #             b[k,j] <- (rgamma(1,shape=(nu[[i]]+2*p)/2,rate=c(t(a[k,,j]-q[j,,i])%*%solve(Lambda[[i]])%*%(a[k,,j]-q[j,,i])+d[j,i]+nu[[i]])/2))
        #           
        #           
        #           e2[j,i] <- mean(b[,j])
        #           e1[j,i] <- mean(log(b[,j]))
        #           
        #           bla1 <- matrix(NA,50,p)
        #           bla2 <- array(NA,c(50,p,p))
        #           
        #           for(k in 1:50)
        #           { 
        #             bla1[k,] <- c(a[k,,j] * b[k,j])  
        #             bla2[k,,] <- b[k,j] * a[k,,j] %*% t(a[k,,j])          
        #           }
        #           
        #           e3[j,,i] <- apply(bla1, 2, mean)
        #           e4[j,,,i] <- apply(bla2, 2:3, mean)
        stop("Kein MCMC")
      }
      
      
    } # for schleife über gruppen
    
    
    #updating parameters
    
    # orig     
    #       mutmp = 0
    #       mutmp1 = 0
    #       for(k in 1:n){
    #         mutmp = mutmp + pro2[k,i]*e2[k,i] * y[k,] 
    #         mutmp1 = mutmp1 + pro2[k,i]*e3[k,,i]
    #       }
    #neu
    mutmp <- colSums(e2[,i] * y)
    mutmp1 <- colSums(e3[,,i])
    
    # so gelassen      
    mu_neu[[i]] = (mutmp - Delta[[i]] %*% mutmp1) / sum(e2[,i])
    
    # orig      
    #       deltmp1 = 0
    #       deltmp2 = 0
    #       for(k in 1:n){
    #         deltmp1 = deltmp1 +  pro2[k,i] * e3[k,,i] %*% t(y[k,] - mu_neu[[i]])
    #         deltmp2 = deltmp2 + pro2[k,i] * e4[k,,,i]
    #       }
    
    # neu
    YneuCenter <- t(t(y) - c(mu_neu[[i]])) 
    deltmp1 <-   t((e3[,,i])) %*% YneuCenter
    
    deltmp2 <-  blockSum(e4list, n=n)
    #Reduce("+", e4list)
    
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
    #       sigtmp = 0
    #       for(k in 1:n){
    #         cenmu <- y[k,] - mu_neu[[i]]
    #         sigtmp  = sigtmp +  pro2[k,i] *(cenmu) %*% t(cenmu)  * e2[k,i] 
    #       }
    
    #neu
    sigtmp <-   t(YneuCenter  * e2[,i]) %*%  (YneuCenter)
    # orig 
    sigtmp = sigtmp + diagDel %*% deltmp2  %*%diagDel -
      diagDel %*% deltmp1 - t(deltmp1) %*% diagDel
    
    ind <- lower.tri(sigtmp)
    sigtmp[ind] <- t(sigtmp)[ind]
    
    
    Sigma_neu[[i]] <- 1 / sum(pro2[,i]) * sigtmp
    
    nu_neu[[i]]<-uniroot(function(nu) nuf(nu, pro2[,i], e2[,i]/pro2[,i], e1[,i]/pro2[,i]),interval=c(0.1,10e7))$root    
    
    iter <- iter+1
    if(verbose) cat("Iteration: ", iter, "\nrelative error", abs(lik_ges / lik_ges_neu-1), "\nLoglikelihood: ", lik_ges_neu, unlist(ttsigma),
                    "\npro", unlist(pro1), "\nmu", unlist(mu_neu), "\nSigma", unlist(Sigma_neu), "\ndelta", unlist(delta_neu), "\nnu", unlist(nu_neu), "\n\n")
    
  } # while schleife zu
  
  
  pro1<-apply(pro2,2,sum)/n
  
  
  
  #lik_ges <- lik_ges_neu
  
  #lik_ges_neu <- sum(log(dmixst(y,list(pro=pro1,mu=mu_neu,Sigma=Sigma_neu,delta=delta_neu,nu=nu_neu))))
  
  
  
  
  
  
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

