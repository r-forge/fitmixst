#' An S4 class the stores the parameters of a mixed skew t distribution
#' @name mixpara-class
#' @aliases mixpara-class
#' @docType class
#' @keywords classes
#' @export
setClass("mixpara",representation(pro="numeric",mu="list",
				Sigma="list",delta="list",
				nu="list", p = "integer", g="integer"),
		validity=function(object){
			if(length(object@pro)==object@g && 
					length(object@mu)==length(object@Sigma) && 
					length(object@Sigma)==length(object@delta) && 
					length(object@delta) ==length(object@nu)) TRUE
			else FALSE})



#' An S4 class the stores the S4 class mixpara, the input data, logLik, S2, S3, iter and call
#' @name fitmixout-class
#' @aliases fitmixout-class
#' @docType class
#' @keywords classes
#' @export
setClass("fitmixout",representation(y="matrix",X = "matrix",beta = "matrix", resid ="matrix", para="mixpara",logLik="numeric",
				iter="numeric",call="call",	empcov="matrix", posteriori = "matrix", likConvergence = "numeric"))

