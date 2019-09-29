#' Generates random observation for the normal distribution
#'
#' @param n number of random observations.
#' @param p dimension of the normal distribution.
#' @param mean a mean vector with length p.
#' @param cov a covariance matrix of dimension pxp.
#' @return gernerates random observations for the normal distribution.
#' @keywords random observations normal distribution.
#' @export
#' @examples
#' 
#' mu=1; Sigma=1;
#' y <- rmvn(n=100,p=1,mean=mu,cov=Sigma)

rmvn<-function(n, p, mean = rep(0, p), cov = diag(p)) 
{
	cov <- as.matrix(cov)
	if (nrow(cov) != ncol(cov)) {
		stop("cov must be a square matrix")
	}
	if (length(mean) != nrow(cov)) {
		stop("mean and cov have non-conforming size")
	}
	rmvnorm(n, mean = mean, sigma = cov, method = "chol")
}

