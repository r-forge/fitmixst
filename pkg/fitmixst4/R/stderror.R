#' The value of the mixed densities of multivarite skew t distributions.
#'
#' @param x grid for which the standard errors are calculated
#' @param model an object of the parameters of the mixture of multivariate skew t distributions.
#'pro a vector for the mixture ratios, mu the values, sigma, delta, nu
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @return numeric vector with the values of the mixture of the density functions.
#' @keywords standard error function
#' @useDynLib fitmixst4
#' @export
#' @examples
#' pro = 1; mu=1; Sigma=3; delta=2; nu=10;
#' para = list(pro=pro,mu=mu,Sigma=Sigma,delta=delta,nu=nu)
#' y=rmixst(100,para)
#' out <- fitmixst(y,g=1,method="kmeans")
#' std = stderror(y, out)


stderror <- function (x, model, log = FALSE){
  if(is(model,"fitmixout")){
    if(model@para@p == 1){
      all = c(model@para@pro, unlist(model@para@mu), unlist(model@para@Sigma), unlist(model@para@delta), 
              unlist(model@para@nu))
      empcov = model@empcov
      gg = model@para@g
    }
    else{
      cat("error: dimension > 1")
      return()
    }
  }
  else
  {
    cat("error: input not a fitmixout object")
  }
  .Call("stderrorC", x, all, empcov, gg, log, PACKAGE = "fitmixst4")
}



