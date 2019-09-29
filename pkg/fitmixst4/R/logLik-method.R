#' extract plot
#' 
#' @name logLik-method
#' @aliases logLik,fitmixout-method
#' @docType methods
#' @importFrom stats4 logLik
#' @exportMethod logLik
setMethod("logLik","fitmixout",
          function(object){
            if(sum(dim(object@beta)) == 0){val <-object@logLik
                           #attr(val, "df") <- 5 * length(object@para@pro)
                           attr(val, "df") <- (2*object@para@p + 1) * object@para@g - 1 + 
                             object@para@g*object@para@p*(object@para@p +1) / 2 + object@para@g
                           attr(val, "nall") <- length(object@y) / object@para@p
                           attr(val, "nobs") <- length(object@y) / object@para@p
                           class(val) <- "logLik"
                           val}
            else{
              val <-object@logLik
              #attr(val, "df") <- 5 * length(object@para@pro)
              #regresseion only 1dim
              attr(val, "df") <- (2*object@para@p + 1) * object@para@g - 1 + 
                object@para@g*object@para@p*(object@para@p +1) / 2 + object@para@g -1+
              #anzahl der betas und 1 para weniger, wegen mu
              nrow(object@beta)
              attr(val, "nall") <- length(object@y) / object@para@p
              attr(val, "nobs") <- length(object@y) / object@para@p
              class(val) <- "logLik"
              val}
                                           })
