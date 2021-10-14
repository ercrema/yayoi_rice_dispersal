cov_GPL2 <- nimbleFunction(
			   run = function(dists = double(2), rhosq = double(0), etasq = double(0), sigmasq = double(0)) {
				   returnType(double(2))
				   n <- dim(dists)[1]
				   result <- matrix(nrow = n, ncol = n, init = FALSE)
				   deltaij <- matrix(nrow = n, ncol = n,init = TRUE)
				   diag(deltaij) <- 1
				   for(i in 1:n)
					   for(j in 1:n)
						   result[i, j] <- etasq*exp(-rhosq*dists[i,j]^2)+sigmasq*deltaij[i,j]
				   return(result)
			   })
