## function to perform confidence interval check given:
## 1. weight function mu
## 2. sigma function for BM
## 3. ??jump function for jump
check.CI <- function(mu, sigma, jump, plotpoints, EBAY, Velo, N.sim, CI){
	## jump is not included yet...
	## perform theoratical value check 
	prop <- matrix(0,length(EBAY), length(CI))
	delta <- 7/(length(plotpoints)-1)
	for(i in 1:length(EBAY)){
	    time <- EBAY[[i]]$Time
	    logprice.r <- log(EBAY[[i]]$Price)
	    
		## save the floor index for each bid time
		## and the coefficient to multiply on floor time
		## coef to multiply on ceiling time is 1-coef_floor
		index <- rep(0, length(time))
		time.floor <- rep(0, length(time))
		coef <- rep(0, length(time))
		for(j in 1:length(logprice.r)){
			# note here +1 is because vector ordering starts from 1
			index[j] <- as.integer(time[j]/delta) + 1
			time.floor[j] <- (index[j] - 1) * 7/(length(plotpoints) - 1)
			coef[j] <- (time[j] - time.floor[j])/delta			
		}
		
		logprice<- matrix(0, N.sim, length(time))
		for(j in 1:N.sim){
			logprice.all <- simApath(mu,sigma,jump=NULL,plotpoints,realpath = EBAY[[i]]$Price, velo.initial = Velo[,1])
			logprice[j, ] <- logprice.all[index] * coef + logprice.all[index+1] * (1-coef)			
		}
		
		# check CI from here
		# now we have real price and sim price
		for(j in 1:length(CI)){
			qt <- c(CI[j]/2, 1- CI[j]/2)
			band.ci <- apply(logprice, 2, function(x){quantile(x, qt)})
			## if within band, tag.ci should be -1 or 0
			tag.ci <- sign((logprice.r - band.ci[1,])*(logprice.r - band.ci[2,]))
			count.out <- length(which(tag.ci == 1))
			prop[i, j] <- 1 - count.out/length(time)
		}
		
		
		
	}
	
	return(prop)
}

simApath <- function(mu, sigma, jump, plotpoints, realpath, velo.initial){
	path <- rep(0, length(plotpoints))
	delta <- 1/(length(plotpoints)-1)
	path[1] <- log(realpath[1])
	path[2] <- path[1] + sample(size = 1, x = velo.initial)
	 
	for(i in 3:length(plotpoints)){
		## Here add the min change being 0??
		## Not sure if leads to some theoretical difficulty
		path[i] <- path[i-1] + (path[i-1] - path[i-2])*(1 + mu[i] * delta + sigma[i] * rnorm(1, 0, sqrt(delta)))
	}
	path
}