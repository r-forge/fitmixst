#' extract probability
#'
#' @name predict-method
#' @aliases predict,fitmixout-method
#' @docType methods
#' @importFrom stats predict
#' @exportMethod predict
setMethod("predict","fitmixout",
          function(object){
            if(sum(dim(object@beta)) == 0){
            dmixst(object@y,
                   list(pro=object@para@pro,mu=object@para@mu,
                        Sigma=object@para@Sigma,delta=object@para@delta,
                        nu=object@para@nu))}
            else{
              object@X %*% object@beta
            }
          })

