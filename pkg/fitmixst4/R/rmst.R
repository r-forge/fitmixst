#' Generates random observation for the skew t distribution
#'
#' @param n number of random observations.
#' @param p dimension of the skew t distribution.
#' @param mean a mean vector with length p.
#' @param cov a covariance matrix of dimension pxp.
#' @param del a skew parameter of the length p.
#' @param nu the degree of freedoms. 
#' @return gernerates random observations for the skew t distribution.
#' @keywords random observations skew t distribution.
#' @export
#' @examples
#' 
#' mu=1; Sigma=1; delta=3; nu=3;
#' y <- rmst(n=100,p=1,mean=mu,cov=Sigma,nu=nu,del=delta)

rmst<-function(n, p, mean = rep(0, p), cov = diag(p), nu = 10, del = rep(0, p)) 
{
	g <- rgamma (n, nu/2, nu/2)
	z = matrix(abs(rnorm(p*n, 0,1) / sqrt(g)),ncol=p)
	x <- t ( t ( rmvn (n, p, cov = cov) / sqrt (g)) + c(mean))
  if(p>1)
    del = diag(c(del))
	as.matrix(z %*% del + x)
}


