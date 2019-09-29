#' Internal fitting function
#'
#' @param y a multidimensional input vector.
#' @param g the number of groups.
#' @param error the desired relative error. (default 1e-5)
#' @param itermax the maximum of iterations. (default 1000)
#' @param pro the number of groups.
#' @param mu the desired relative error. (default 1e-5)
#' @param Sigma the maximum of iterations. (default 1000)
#' @param delta of finding initial values. (default "kmeans")
#' @param nu of finding initial values. (default "kmeans")
#' @param verbose of finding initial values. (default "kmeans")
#' @return a object of the class fitmixst.
#' @keywords fit mixed skew t.
#' @export



mmixstEstC <- function (y, g, itermax, error, pro, mu, Sigma, delta, nu, verbose=F) 
  .Call("fitmixstMVC", y, g, itermax, error, pro, mu, Sigma, 
        delta, nu, verbose, PACKAGE = "fitmixst4")
