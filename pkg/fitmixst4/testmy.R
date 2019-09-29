library(fitmixst4)
library(mnormt)
require(rbenchmark)
n <- 100
y <- cbind(rnorm(n), rnorm(n,1))

set.seed(1);
system.time(
out1 <- fitmixst(y,g=1,method="kmeans",itermax=200, rel.error=0)
)
set.seed(1);
system.time(
  out2 <- fitmixst(y,g=1,method="kmeans", useR = T,itermax=200, rel.error=0, R2=FALSE,
                   verbose=F)
)
fix(fitmixst) 
system.time(
  out3 <- fitmixst(y,g=1,method="kmeans", useR = T, R2=TRUE, itermax=200, rel.error=0,
                   verbose=F)
)

benchmark(replications=10,
          fitmixst(y,g=1,method="kmeans", useR = T,itermax=200, rel.error=1e-5, R2=FALSE,
                   verbose=F),
          fitmixst(y,g=1,method="kmeans", useR = T, R2=TRUE, itermax=200, rel.error=1e-5,
                   verbose=F),
          columns=c('test', 'elapsed', 'relative'))

Rprof("bla.out")
#fitmixst(y,g=1,method="kmeans", useR = T,itermax=20, rel.error=1e-5, andi=FALSE,
#                 verbose=F)
fitmixst(y,g=1,method="kmeans", useR = T, R2=TRUE, itermax=20, rel.error=1e-5,
         verbose=F)
Rprof(NULL)
summaryRprof("bla.out")

ttsigma <- rep(list(diag(3),diag(3)*2,diag(3)*3),100)
ttmu2 <- rep(list(1:3),100)
ttmu <- 1:(3*100)

ttmulist <- lapply(seq_len(nrow(ttmu)), function(ii) ttmu[ii,])





bla <- function()
{
scal <- matrix(rep(ttmu, each=9), ncol=3, byrow=T)
M2 <- do.call(rbind, ttsigma)
               
M2* scal
}
mapply("*", ttsigma, as.list(ttmu), SIMPLIFY=F)
ttmu <- 1:(100*100)
benchmark(replications=1e4,
          t(sapply(ttmu2, function(x) x)),
          matrix(unlist(ttmu2), ncol=2, byrow=T),
          do.call(rbind, ttmu2),
          columns=c('test', 'elapsed', 'relative'))

system.time(pmt.vec(mat=M2, mean=c(-1,.1,.1), S=diag(3),df=4))
system.time(pmt.vec2(mat=M2, mean=c(-1,.1,.1), S=diag(3),df=4))

M2 <- M2[1:3,]
benchmark(replications=10,
          #pmt.vec(mat=M2, mean=c(-1,.1,.1), S=diag(3),df=4),
          pmt.vec2(mat=M2, mean=c(-1,.1,.1), S=diag(3),df=4),
          pmt.vec3(mat=M2, mean=c(-1,.1,.1), S=diag(3),df=4),
          #pmt.vec4(mat=M2, mean=c(-1,.1,.1), S=diag(3),df=4),
          #pmt.vec5(mat=M2, mean=c(-1,.1,.1), S=diag(3),df=4),
          columns=c('test', 'elapsed', 'relative'))

y <- rbind(rmvnorm(100, mean=c(0,0)), rmvnorm(100, mean=c(2,2)))
xg <- seq(min(y[,1]), max(y[,1]), length=100)
yg <- seq(min(y[,2]), max(y[,2]), length=100)
zgrid <- expand.grid(xg,yg)
dmix <- function(x) .5 * dmvnorm(x) + .5 *dmvnorm(x, c(2,2))
z <- matrix(dmix(zgrid), 100, byrow=T)
plot(y)
contour(xg,yg, z, add=T)
