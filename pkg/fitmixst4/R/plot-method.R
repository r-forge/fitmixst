#' extract plot
#' 
#' @name plot-method
#' @aliases plot,mixpara-method
#' @docType methods
#' @rdname plot-method
#' @importFrom graphics plot
#' @exportMethod plot
setMethod("plot","mixpara",function(x,y,xlab="x-axis",ylab="density",from=-1,to=1,...){
  if(missing(y)){
    if(x@p == 1){
    curve(dmixst(z,x),xname="z",from=from,to=to,xlab=xlab,ylab=ylab,...)}
  }
  else cat("Error; use arguments from and to in plot function")
})



#' @aliases plot,fitmixout-method
#' @docType methods
#' @rdname plot-method
#' @importFrom graphics plot
#' @exportMethod plot
setMethod("plot","fitmixout",function(x,y,xlab="x-axis",ylab="density",main = x@call, type="l", from =min(x@y), to = max(x@y), ...){
  if(missing(y)){
    if(x@para@p == 1){
      curve(dmixst(z,x@para),xname="z",from=from,to=to,type=type , xlab=xlab,ylab=ylab, main = main, ...)
    }
    
    if(x@para@p == 2){
      n = 100
    x1 <- seq(min(x@y[,1]) , max(x@y[,1]), length.out=n)
    x2 <- seq(min(x@y[,2]) , max(x@y[,2]), length.out=n)
    
    x1x2 <- expand.grid(x1,x2)
    
    z <- matrix(dmixst(x1x2, x@para), ncol = n, byrow=F)
    
    
#     z <- matrix(0, n, n)
#     for(i in 1:n){
#       for(j in 1:n){
#         z[i,j] <- dmixst(c(x1[i], x2[j]), x@para)
#       }
#     }
    contour(x1,x2,z,...)
    
  }
  }
  else cat("Error; use arguments from and to in plot function")
})


#' @aliases lines,fitmixout-method
#' @docType methods
#' @rdname plot-method
#' @importFrom graphics lines
#' @exportMethod lines
setMethod("lines","fitmixout",function(x,xlab="x-axis",ylab="density",...){
  if(x@para@p == 1){
  lines(sort(x@y),dmixst(sort(x@y),x@para),...)
  }})

#' @aliases points,fitmixout-method
#' @docType methods
#' @rdname plot-method
#' @importFrom graphics points
#' @exportMethod points
setMethod("points","fitmixout",function(x,...){
  if(x@para@p == 1){
  points(sort(x@y),dmixst(sort(x@y),x@para),...)}
})







