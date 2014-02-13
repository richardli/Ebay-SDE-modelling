## combine all the residuals after taking out determinant part
require(rootSolve)
library(rootSolve)
FitsubSDE.jump <- function(index){
	all.res <- NULL
	Acc.temp <- Acc.smooth[index, ]
	Velo.temp <- Velo.smooth[index, ]
	mu.temp <- fitSDE_mu(index, kk)
	for(i in 1:length(plotpoints)){
		temp <- (Acc.temp[, i] - mu.temp[i]*Velo.temp[, i])/Velo.temp[,i]
		all.res <- c(all.res,remove_outliers(temp))
	}	
	all.res <- all.res[-which(is.na(all.res))]
	## calculate Cumulants of the residuals
	m1 <- moment(all.res, order = 1)
	m2 <- moment(all.res, order = 2)
	m3 <- moment(all.res, order = 3)
	m4 <- moment(all.res, order = 4)
	k1 <- m1
	k2 <- m2 - m1 ^ 2
	k3 <- m3 - 3 * m1 * m2 + 2 * m1 ^ 3
	k4 <- m4 - 4*m1*m3 - 3*m2^2 + 12*m2*m1^2 - 6*m1^4
	k4 <- all.cumulants(all.moments(x, order.max = 4))[5]
	
	
	func <- function(x) x^4 - 2* k3/k1 *x^2 + 3*k4/2/k1* x - k3/2/k1
	mu <- uniroot.all(func, c(-10,100))
	mu <- mu[which(mu > 0)]
	lambda <- k1 / mu
	ita2 <- (k3 - mu^2 *k1)/(3*k1)
	sigma2 <- k2 - k1 * mu  + ita2
	
	return(c(mu, lambda, ita2, sigma2))
}

