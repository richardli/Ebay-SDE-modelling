## function to perform confidence interval check given:
## 1. weight function mu
## 2. sigma function for BM
## 3. ??jump function for jump
check.CI <- function(mu, sigma, jump, plotpoints, EBAY, Velo, N.sim ){
	## jump is not included yet...
	## perform theoratical value check 
	for(i in 1:length(EBAY)){
		for(j in 1:N.sim){
			logprice <- simApath(mu,sigma,jump=NULL,plotpoints,realpath = EBAY[[i]]$Price, velo.initial = Velo[i,1])
			
			
			## insert codes to calculate CI here
		
		}
	}
	
	
}

simApath <- function(mu, sigma, jump, plotpoints, realpath, velo.initial){
	path <- rep(0, length(plotpoints))
	delta <- 1/(length(plotpoints)-1)
	path[1] <- log(realpath[1])
	path[2] <- path[1] + sample(size = 1, x = velo.initial) #velo.initial * delta 
	for(i in 3:length(plotpoints)){
		path[i] <- path[i-1] + (path[i-1] - path[i-2])*(1 + mu[i] * delta + sigma[i] * rnorm(1, 0, sqrt(delta)))
	}
	path
}