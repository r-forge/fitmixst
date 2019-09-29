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
#' pro = 1; mu=1; Sigma=3; delta=2; nu=3;
#' para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
#' y <- rmixst(100,para)
#' out <- fitmixst(y,g=1,method="kmeans")

regSkewt <- function(y, X, rel.error = 1e-3, itermax=1000, 
                     method = "kmeans", a = 0.5, para = NULL,verbose=F, mcmc=F, 
                     useR = F, R2 = FALSE, ...)
{
  
  if(mcmc)
    warning("load the package tmvtnorm manually!")
  
  d1 <- dim(y)

  if(is.null(d1) || sum(dim(d1)) < 2){
    
    y <- as.vector(y)
    g <- 1
    
  beta <- vector(length=ncol(X))
    
    
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
      }
      else stop("method and input parameters are not consistent")
      
    }
    
    tryCatch(est <- regtest(y=y, X=X, g=g, itermax=itermax, error=rel.error, 
                             pro=pro, mu=mu, Sigma=Sigma, delta=delta, 
                             nu=nu, verbose = ifelse(verbose,1,0)), "C++Error" = function(e){
                               cat( sprintf( "C++ exception of class '%s' : %s\n", class(e)[1L], e$message  ) )
                             } )
    
    
  }
  else{
      
      
    if(method =="kmeans"){
      
      #a <- 0.5
      g <- 1
      init <- kmeans(y, g,nstart=5)
      
      pro <- init$size/length(y[,1])
      mu <- delta <- Sigma <- nu <- vector(mode="list", length=g)
      
      for (j in 1:g) {
        if(sum(init$cluster == j) > 1){
          Sigma[[j]] <- var(y[init$cluster == j,]) + (a-1)*diag(var(y[init$cluster == j,]))
          mu[[j]] <- init$centers[j, ]
          delta[[j]] <- sign(sum((y[init$cluster == j,] -(rep(mu[[j]], sum(init$cluster == j))))^3)) * 
            sqrt((1-a)* diag(var(y[init$cluster == j,]))/(1-2/pi))
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
          warning("Parameter g to large. Too many groups specified!")
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
      
      tryCatch(est <- regtestMV(y=y, X=X, g=g, itermax=itermax, error=rel.error, 
                               pro=pro, mu=mu, Sigma=Sigma, delta=delta, 
                               nu=unlist(nu), verbose = ifelse(verbose,1,0)), "C++Error" = function(e){
                                 cat( sprintf( "C++ exception of class '%s' : %s\n", class(e)[1L], e$message  ) )
                               } )
    
    
  }
    
    
    #print(est$mu)
    est$mu <- list(est$mu)
    est$Sigma <- list(est$Sigma)
    est$delta <- list(est$delta)
    est$nu <- list(est$nu)
    est$posteriori <- matrix(unlist(est$posteriori),ncol=1)
    
  
  est$call <- sys.call()
  
  para <- new("mixpara",pro=est$pro,mu=est$mu,Sigma=est$Sigma,
              delta=est$delta,nu=est$nu, p=est$p, g=as.integer(est$g))
  out <- new("fitmixout", y=as.matrix(y), X=X, beta = est$beta, resid = y - X%*%est$beta, para=para, logLik=est$logLik, 
             iter=est$iter,call=est$call, empcov=est$empcov, posteriori = est$posteriori)
  out
  
}
