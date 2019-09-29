#' Fitted mixed model
#' @usage fitmixst(y, g, rel.error = 1e-3, itermax=1000, method = "kmeans", para = NULL, verbose = F, mcmc = F, useR = F, ...)
#' @param y a multidimensional input vector.
#' @param g the number of groups.
#' @param rel.error the desired relative error. (default 1e-3)
#' @param itermax the maximum of iterations. (default 1000)
#' @param method of finding initial values. (default "kmeans")
#' @param para the initial input parameters. Note: if you choose input parameters the kmeans method is disabled.
#' para a list of the input parameters. para = list (pro, mu, Sigma, delta, nu).
#' @param verbose debug mode.
#' @param mcmc calculates moments moments with mcmc in R
#' @param useR use algorithm programmed in R (default is C++ Algorithm).
#' @param R2 use alternative R code.
#' @param ... other paramters.
#' @useDynLib fitmixst4
#' @export
#' @return an S4 object of the class fitmixout.
#' @keywords fit skew t distribution.
#' @examples
#' 
#' pro = 1; mu=1; Sigma=3; delta=2; nu=10;
#' para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
#' y <- rmixst(100,para)
#' out <- fitmixst(y,g=1,method="kmeans")

fitmixst <- function(y, g, rel.error = 1e-3, itermax=1000, 
                     method = "kmeans", a = 0.5, para = NULL,verbose=F, mcmc=F, 
                     useR = F, R2 = FALSE, beta = 1,...)
{
  
  if(!(method %in% c("kmeans","self")))
    stop("Method should be either kmeans or self.")
  
  
  if(mcmc)
    warning("Load the package tmvtnorm manually!")
  
  d1 <- dim(y)
  
  if(!is.null(d1)){
    
    if(method =="kmeans"){
      
      #a <- 0.5
      init <- kmeans(y, g,nstart=5)
      
      pro <- init$size/length(y[,1])
      mu <- delta <- Sigma <- nu <- vector(mode="list", length=g)
      
      for (j in 1:g) {
        if(sum(init$cluster == j) > 1){
          Sigma[[j]] <- cov(y[init$cluster == j,]) + (a-1)*diag(diag(cov(y[init$cluster == j,])))
          
          while(min(eigen(Sigma[[j]])$values) <= 0 && a < 0.9){
            a = round(a+0.1,2)
            Sigma[[j]] <- cov(y[init$cluster == j,]) + (a-1)*diag(diag(cov(y[init$cluster == j,])))
            warning("Parameter a changed for positive definite starting value of sigma.")
          }
          
          if(min(eigen(Sigma[[j]])$values) <= 0){
            Sigma[[j]] <- cov(y[init$cluster == j,])
            warning("Scale matrix sigma changed for positive definite starting value.")
          }
          
          if(min(eigen(Sigma[[j]])$values) <= 0){
            Sigma[[j]] <- diag(diag(Sigma[[j]]))
            warning("Scale matrix sigma changed for positive definite starting value.")
          }
          
          mu[[j]] <- init$centers[j, ]
          delta[[j]] <- sign(sum((y[init$cluster == j,] -(rep(mu[[j]], sum(init$cluster == j))))^3)) * 
            sqrt((1-a)* diag(cov(y[init$cluster == j,]))/(1-2/pi))
          mu[[j]] <- init$centers[j, ] - sqrt(2/pi) * delta[[j]]
          dimnames(Sigma[j]) <- NULL
          names(mu[j]) <- NULL
          names(delta[j]) <- NULL
          nu[[j]] <- 40
        }else{
          Sigma[[j]] <- diag(d1[2])
          mu[[j]] <- init$centers[j, ]
          delta[[j]] <- matrix(0, nrow=d1[2], ncol =1)
          mu[[j]] <- init$centers[j, ]
          dimnames(Sigma[j]) <- NULL
          names(mu[j]) <- NULL
          names(delta[j]) <- NULL
          nu[[j]] <- 40
          warning("Parameter g to large. No observations in one group with kmeans initialization. Too many groups specified!")
        }
      }	
      
    }
    if(method =="self" && ! is.null(para)){
      if(is(para,"mixpara")){
        pro  <- para@pro
        mu <- para@mu
        Sigma  <- para@Sigma
        delta  <- para@delta
        nu <- para@nu
        p <- para@p
        g <- para@g
      }
      else if(is(para,"fitmixout")){
        pro  <- para@para@pro
        mu <- para@para@mu
        Sigma  <- para@para@Sigma
        delta  <- para@para@delta
        nu <- para@para@nu
        p <- para@para@p
        g <- para@para@g
      }
      else
      {
        pro<- para$pro
        mu <- para$mu
        Sigma <- para$Sigma
        delta <- para$delta
        nu <- para$nu
        g <- length(para$pro)
        if(g == 1){
          p <- length(mu[[1]])
        }
        if(g > 1){
          p <- length(mu[[1]])
        }
      }
    }
    
    if(mcmc || useR){
      if(R2)
      {
        tryCatch(est <- mmixstEst2(y, g, itermax = itermax, error = rel.error,pro=pro, mu=mu, Sigma=Sigma, delta=delta, 
                                   nu=nu, verbose = verbose, mcmc=mcmc), "mmixstEst error" = function(e){
                                     cat( "Error in mmixstEst" )}) 
      }
      if(!R2)
      {
        tryCatch(est <- mmixstEst(y, g, itermax = itermax, error = rel.error,pro=pro, mu=mu, Sigma=Sigma, delta=delta, 
                                  nu=nu, verbose = verbose, mcmc=mcmc), "mmixstEst error" = function(e){
                                    cat( "Error in mmixstEst" )})  
      }
      
    }else{
      
      tryCatch(est <- mmixstEstC(y=y, g=g, itermax=itermax, error=rel.error, 
                                 pro=pro, mu=mu, Sigma=Sigma, delta=delta, 
                                 nu=unlist(nu), verbose = ifelse(verbose,1,0)), "C++Error" = function(e){
                                   cat( sprintf( "C++ exception of class '%s' : %s\n", class(e)[1L], e$message  ) )
                                 } )
      
    } 
  }else{
    
    
    y <- as.vector(y)
    g <- as.integer(g)
    
    
    
    #using kmeans to get inital values
    if(method == "kmeans" || is.null(para)){
      
      #a  <- 0.9
      init <- kmeans(y, g,nstart=5)
      
      pro <- init$size/length(y)
      mu <- delta <- Sigma <- vector(length=g)
      for (j in 1:g) {
        Sigma[j] <- a*var(y[init$cluster == j])
        delta[j] <- sign(sum((y[init$cluster == j] -(init$centers[j, ]))^3)) * 
          sqrt((1-a)*var(y[init$cluster == j])/(1-2/pi))
        mu[j] <- init$centers[j, ] - sqrt(2/pi) * delta[j]
        dimnames(Sigma[j]) <- NULL
        names(mu[j]) <- NULL
        names(delta[j]) <- NULL
      }	
      nu <-rep(40,g)
      
      
    }
    else{
      
      if(mcmc)
        stop("mcmc only possible for p > 1")
      
      if(method =="self" && ! is.null(para)){
        if(is(para,"mixpara")){
          pro  <- para@pro
          mu <- para@mu
          Sigma  <- para@Sigma
          delta  <- para@delta
          nu <- para@nu
          p <- para@p
          g <- para@g
        }
        else if(is(para,"fitmixout")){
          pro  <- para@para@pro
          mu <- para@para@mu
          Sigma  <- para@para@Sigma
          delta  <- para@para@delta
          nu <- para@para@nu
          p <- para@para@p
          g <- para@para@g
        }
        else
        {
          pro<- para$pro
          mu <- para$mu
          Sigma <- para$Sigma
          delta <- para$delta
          nu <- para$nu
          g <- length(para$pro)
          if(g == 1){
            p <- length(mu[[1]])
          }
          if(g > 1){
            p <- length(mu[[1]])
          }
        }
        
        mu = unlist(mu); Sigma = unlist(Sigma); delta = unlist(delta); nu = unlist(nu)
      }
      else stop("method and input parameters are not consistent")
      
    }
    
    tryCatch(est <- mixstEst(y=y, g=g, itermax=itermax, error=rel.error, 
                             pro=pro, mu=mu, Sigma=Sigma, delta=delta, 
                             nu=nu, verbose = ifelse(verbose,1,0), beta = beta), "C++Error" = function(e){
                               cat( sprintf( "C++ exception of class '%s' : %s\n", class(e)[1L], e$message  ) )
                             } )
    
    
    est$mu <- as.list(est$mu)
    est$Sigma <- as.list(est$Sigma)
    est$delta <- as.list(est$delta)
    est$nu <- as.list(est$nu)
    est$posteriori <- matrix(unlist(est$posteriori),ncol=1)
    
    
  }
  est$call <- sys.call()
  
  #sort by size
  
  if(g > 1){
    orderOut <- order(sapply(est$mu, "[[", 1))
    est$mu <- est$mu[orderOut]
    est$Sigma <- est$Sigma[orderOut]
    est$delta <- est$delta[orderOut]
    est$nu <- est$nu[orderOut]
    est$pro <- est$pro[orderOut]
    numParam = 2*est$p + est$p*(est$p +1) / 2 + 1
    numPi = est$g - 1
    
    numPro <- numMu <- numSigma <- numDelta <- numNu <- vector(mode="list",length=est$g)
    numMu = list()
    numPro = seq(1, est$g, length.out=est$g-1)
    for(k in 1:est$g){
      numMu[[k]] = est$g-1 + 1:est$p + est$p * (k-1)
      numSigma[[k]] = est$g-1 + est$p * est$g + 1:(est$p*(est$p +1) / 2) + est$p*(est$p +1) / 2 *(k-1)
      numDelta[[k]] = est$g-1 + est$p * est$g +est$p*(est$p +1) / 2 * est$g + 1:est$p + est$p * (k-1)
      numNu[[k]] = est$g-1 + est$p * est$g +est$p*(est$p +1) / 2 * est$g + est$p * est$g + k
    }
    
    #order numPro not correct
    orderOutMatrix = c(numPro, unlist(numMu[orderOut]), unlist(numSigma[orderOut]), 
                       unlist(numDelta[orderOut]), unlist(numNu[orderOut]))
    
    
    est$empcov <- est$empcov[orderOutMatrix,orderOutMatrix]
  }
  
  if(g>1)
    est$pro[1] = 1 - sum(est$pro[-1])
  
  para <- new("mixpara",pro=est$pro,mu=est$mu,Sigma=est$Sigma,
              delta=est$delta,nu=est$nu, p=est$p, g=as.integer(est$g))
  out <- new("fitmixout", y=as.matrix(y), para=para, logLik=est$logLik, 
             iter=est$iter,call=est$call, empcov=est$empcov, posteriori = est$posteriori, likConvergence = est$likConvergence)
  
  out
  
}
