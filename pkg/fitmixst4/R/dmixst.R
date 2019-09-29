#' The value of the mixed densities of multivarite skew t distributions.
#'
#' @param x a vector.
#' @param para an object of the parameters of the mixture of multivariate skew t distributions.
#'pro a vector for the mixture ratios, mu the values, sigma, delta, nu
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return numeric vector with the values of the mixture of the density functions.
#' @keywords density function
#' @useDynLib fitmixst4
#' @export
#' @examples
#' x=c(1,2)
#' pro = 1; mu=1; Sigma=3; delta=2; nu=10;
#' para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
#' dmixst(x,para)


dmixst <- function (x, para, log = FALSE){
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
  
  #univariate
	if(p == 1){
	.Call("dmixstC", x, pro, unlist(mu), unlist(Sigma), unlist(delta), unlist(nu), log, PACKAGE = "fitmixst4")
  
  #multivariate
  }else{
    
    #multivariate density function

      
#       pro <- para$pro
#       mu <- para$mu
#       Sigma <- para$Sigma
#       delta <- para$delta
#       nu <- para$nu
      
      p <- length(mu[[1]])
      #np <- dim(data.frame(x))
      n <- nrow(x)#length(x) / p
      #p <- np[2]
      x <- as.matrix(x)
      x <- matrix(x, n,p)
      
      g <- length(pro)
#       phat <- vector(length = n)
#       
#       for (i in 1:n){
#         phat [i] <- 0
#         
#         for (j in 1:g){ 
#           Omega <- Sigma[[j]] + diag(c(delta[[j]]))^2
#           Omegainv <- solve(Omega)
#           Lambda <- diag(p)-t(diag(c(delta[[j]]))) %*% Omegainv %*% diag(c(delta[[j]]))
#           #cat(t(x[i,] - mu[[j]]))
#           #cat("\n")
#           term1 <- c(sqrt((nu[[j]] + p) / (nu[[j]] + t(x[i,] - mu[[j]]) %*% Omegainv %*% (x[i,] - mu[[j]]))))
#           #cat(term1)
#           mutmp <- diag(c(delta[[j]])) %*% Omegainv %*% (x[i,] - mu[[j]]) * term1
#           
#           phat[i] = phat[i] + pro[i] *  2^p * dmt(x[i,], mean= c(mu[[j]]), S = Omega, df = round(nu[[j]])) * 
#             pmt(x=mutmp, mean= rep(0,p), S = Lambda, df = round(nu[[j]]) + p)
#         }
#       }
#       phat

      .Call("dmmixstC", x, pro, mu, Sigma, delta, unlist(nu), log, PACKAGE = "fitmixst4")
      
      
    
  }
  
  
}


