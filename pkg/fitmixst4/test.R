# TODO: Add comment
# 
# Author: jh
###############################################################################



Sigma <- list(matrix(c(1,2,3,4),2,2))
mu <- list(c(2,2))
delta <- list(c(3,2))
nu <- list(4)
y <- cbind(c(1,2,3,4,5,6),c(6,5,4,3,2,1))
i <- 1

fun <- function(x,y){
	z=5*x;
	j=y/3;
	z+j;
}

myfun <- function(x) mean(rnorm(10000))
system.time(
		bla <- mclapply(1:1e3, myfun, mc.cores=4)
)



require(doSNOW)
registerDoSNOW(makeCluster(4, type = "SOCK"))
getDoParWorkers()
r=100#100replications

MLE<-function()
{
	x=rnorm(1e6)
	y=rnorm(1e6)
	Zeit=system.time(lm(y~x+I(x^2)))[3]
	Zeit=Zeit +system.time(lm(y~x+I(x^2)))[3]
	Zeit=Zeit +system.time(lm(y~x+I(x^2)))[3]
	Zeit=Zeit +system.time(lm(y~x+I(x^2)))[3]
	return(Zeit)
}


out=foreach(i=1:4,.errorhandling="remove", .combine=c) %do% { mle=MLE()}

out2=foreach(i=1:4,.errorhandling="remove", .combine=c) %dopar% { mle=MLE()}

sum(out)/sum(out2)





if(.Platform$OS.type == "windows") {
	detectCores <- function(all.tests = FALSE)
		.Call("ncpus", FALSE, package = "parallel")
} else {
	detectCores <- function(all.tests = FALSE)
	{
		systems <-
				list(darwin  = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null",
						freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
						linux   = "grep processor /proc/cpuinfo 2>/dev/null | wc -l",
						irix    = c("hinv | grep Processors | sed 's: .*::'",
								"hinv | grep '^Processor '| wc -l"),
#             solaris = "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l"
						solaris = "/usr/sbin/psrinfo -p")
		for (i in seq(systems))
			if(all.tests ||
					length(grep(paste("^", names(systems)[i], sep=''), R.version$os)))
				for (cmd in systems[i]) {
					a <- gsub("^ +","", system(cmd, TRUE)[1])
					if (length(grep("^[1-9]", a))) return(as.integer(a))
				}
		NA_integer_
	}
}
