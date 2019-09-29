#' extract summary
#'
#' @name summary-method
#' @aliases summary,fitmixout-method
#' @docType methods
#' @exportMethod summary
setMethod("summary","fitmixout",
          function(object){cat("Fitmixst \n Number of iteration: "); 
                           cat(as.character(object@iter));
                           cat("\n Log Likelhood:"); cat(as.character(round(object@logLik,4)));
                           cat("\n \n"); print(object@para)})
