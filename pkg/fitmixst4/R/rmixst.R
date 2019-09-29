#' Gernates random observation for a mixed skew t distribution
#' @param n number of generated observations.
#' @param para the initial input parameters. 
#' para a list of the input parameters. para = list (pro, mu, Sigma, delta, nu).
#' @return a vector of random observations.
#' @keywords generate random variates.
#' @export
#' @examples
#' 
#' pro = 1; mu=1; Sigma=3; delta=2; nu=3;
#' para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
#' y <- rmixst(100,para)


rmixst<-function(n, para)
{
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
  
  if ( n < g)
  {
    stop("The number of random samples should >= number of components")
  }
   
  #check if mu is a list
  # if mu is not a list => univariate
  if(p==1){
	y <- vector()
	
	if ( abs(sum(pro) -  1) > .Machine$double.eps)
	{
		stop("Sum of pro vector hast to be one")
	}
  pronew = round(n*pro)
  if(g > 1)
    pronew[1] = n - sum(pronew[-1])
	for (i in 1:g)
	{
		y <- c(y,rmst(pronew[[i]],p=1,mean=unlist(mu)[[i]],cov=unlist(Sigma)[[i]],nu=unlist(nu)[[i]],del=unlist(delta)[[i]]))
	}	
  #multivariate generator
  }else{
    p <- length(mu[[1]])
    y <- NULL
    
    if (abs(sum(pro) -  1) > .Machine$double.eps) {
      stop("Sum of pro vector hast to be one")
    }
    pronew = round(n*pro)
    if(g > 1)
      pronew[1] = n - sum(pronew[-1])
    for (i in 1:g) {
      y <- rbind(y, rmst(pronew[[i]], p = p, mean = mu[[i]], cov = Sigma[[i]], 
                         nu = nu[[i]], del = delta[[i]]))
    }
  }
  
  
  
  
  
	return(y)
}

