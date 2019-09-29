
setAs("mixpara", "vector", function(from)
		{
			value <- as.vector(c(from@pro,from@mu,from@Sigma,from@delta,from@nu))
			p <- length(from@pro)
			names(value) = c(rep("pro",p),rep("mu",p),rep("Sigma",p),rep("delta",p),rep("nu",p))
			value
		})
