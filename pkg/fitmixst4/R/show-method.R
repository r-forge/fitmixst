#'extract show
#'
#' @name show-method
#' @aliases show,mixpara-method
#' @docType methods
#' @rdname show-method
#' @importFrom methods show
#' @exportMethod show
setMethod("show","mixpara",function(object){cat("Mixture parameter: \n")
                                            cat("pro \n");print(object@pro); 
                                            cat("mu \n"); print(object@mu);
                                            cat("Sigma \n"); print(object@Sigma);
                                            cat("delta \n"); print(object@delta);
                                            cat("nu \n"); print(object@nu);})



#' @aliases show,fitmixout-method
#' @docType methods
#' @rdname show-method
#' @importFrom methods show
#' @exportMethod show
setMethod("show","fitmixout",function(object){
  if(sum(dim(object@beta)) == 0){
  cat("Fitmixst \n"); print(object@call);
                                              cat("Number of iteration: "); 
                                              cat(as.character(object@iter));
                                              cat("\n Log Likelhood: "); cat(as.character(round(object@logLik,4)));
                                              cat("\n \n"); print(object@para)}
  else{
    cat("Skew t Regression \n"); print(object@call);
  cat("Number of iteration: "); 
  cat(as.character(object@iter));
  cat("\n Log Likelhood: "); cat(as.character(round(object@logLik,4)));
  cat("\n\n beta estimates: "); print(object@beta)
  #cat("\n\nError Distribution");print(object@para)
  }
  })
