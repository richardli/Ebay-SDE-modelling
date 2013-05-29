setwd("~/Dropbox/Ebay_2012/SDE approach")
require("fda")
library("fda")
load("ebay_data.rda") # With elemantary functions
N.auc <- length(EBAY)
##############################################################################################
SPLINE=function(x, y, plotpoints){
	knots <- c(0,1,2,3,4,5,6,6.25,6.5,6.75,6.8125,6.875,6.9375,7)
	numk<-length(knots)
	numplot <- length(plotpoints)
	norder <- 5
	nbasis <- numk + 2
	lambda <- 0.1
	numauct<-N.auc
	wbasis <- create.bspline.basis(rangeval=c(0,7),nbasis=nbasis, norder=5)
	Wfd <- list(0); length(Wfd) <- numauct
	bids=rep(0,length(knots))
	for(j in 1:length(knots)) bids[j]=StepCheck(x,y,knots[j])
	Wfd <- Data2fd(knots,bids,basisobj=wbasis)

	growfdPar <- fdPar(Wfd, 3, lambda)	
	mss <-smooth.basis(knots,bids,growfdPar,wtvec=rep(1,length(knots)))
	xfd <- mss$fd
	splineF <- NULL
	splinef <- NULL
	splineF <- eval.fd(plotpoints,xfd)
	splinef <- eval.fd(plotpoints,xfd,Lfd = 1)
	splineff <- eval.fd(plotpoints, xfd, Lfd = 2)
	return(cbind(splineF,splinef,splineff))
}
Nt <- 2000
plotpoints <- seq(0,7,len = Nt)  #nodes where the smoothed curve is evaluated
LogPrice.smooth <-  matrix(0, N.auc, length(plotpoints))
Velo.smooth <-  matrix(0, N.auc, length(plotpoints))
Acc.smooth <-  matrix(0, N.auc, length(plotpoints))

for(i in 1 : N.auc){
  temp <- SPLINE(EBAY[[i]]$Time,log(EBAY[[i]]$Price), plotpoints)
  LogPrice.smooth[i, ] <- temp[, 1]
  Velo.smooth[i, ] <- temp[, 2]
  Acc.smooth[i, ] <- temp[, 3]
 }
 
par(mfrow = c(2, 2))

plot(plotpoints, LogPrice.smooth[1, ], type = "l", col = "white", ylim = c(-20, 20))
for(k in 1 : N.auc){
	points(plotpoints, LogPrice.smooth[k, ], type = "l", col = k)
}
plot(plotpoints, Velo.smooth[1, ], type = "l", col = "white", ylim = c(-20,10))
for(k in 1 : N.auc){
	points(plotpoints, Velo.smooth[k, ], type = "l", col = k)
}
plot(plotpoints, Acc.smooth[1, ], type = "l", col = "white", ylim = c(-6,6))
for(k in 1 : N.auc){
	points(plotpoints, Acc.smooth[k, ], type = "l", col = k)
}


Acc.M <- apply(Acc.smooth, 2, mean)
Velo.M <- apply(Velo.smooth, 2, mean)
plot(Velo.M, Acc.M, type = "l", xlim = c(0, 1), ylim = c(-.5, .5))
## Here might need to change the index if some auctions are deleted!!
Acc.M.1 <- apply(Acc.smooth[1:93, ], 2, mean)  # Xbox
Velo.M.1 <- apply(Velo.smooth[1:93, ], 2, mean)
Acc.M.2 <- apply(Acc.smooth[94:190, ], 2, mean)  # HP book
Velo.M.2 <- apply(Velo.smooth[94:190, ], 2, mean)
points(Velo.M.1, Acc.M.1, type = "l", lty = 2, col = "red")
points(Velo.M.2, Acc.M.2, type = "l", lty = 3, col = "blue")
######################################################################################
## Simulation
## Linear Trend
par(mfrow=c(1, 2))
Nx <- 2000
delta <- 7/Nx
xx <- seq(0,7,len=Nx)
yy <- xx * 0.4 -1.6
fbasis <- create.fourier.basis(c(0, 7), nbasis = 5)
ftrans <- smooth.basis(xx, yy, fbasis)
plot(xx, eval.fd(xx, ftrans$fd), type = "l")
points(xx, yy, type = "l", lty = 2, col = "red")
## Weight Function
w <- eval.fd(xx, ftrans$fd)
## Simulation
Nsim <- 200
y.sim <- matrix(0, Nsim, Nx)
for(j in 1:Nsim){
	w.noise <- w + rnorm(Nx, mean = 0, sd = 0.01)
	y.sim[j, 1] <- runif(1, -5, 5)
	c <- runif(1, 0, 10) 
	y.sim[j, 2] <- y.sim[j, 1] + abs(rnorm(1, mean = 0, sd = 0.01))
	for (i in 3 : Nx){
	 	y.sim[j, i] <- (w.noise[i] * delta * y.sim[j, i-1] - 2 * y.sim[j, i-1] + y.sim[j, i-2]) / (w.noise[i] * delta -1)
    }
}

plot(xx, y.sim[1, ], type = "l", ylim = c(-5, 15))
for(i in 2 : Nsim){
	points(xx, y.sim[i, ], type = "l", lty = i, col = i)
}
############################################################################################
## Simulation with Brownian Motion
Nsim <- 100
par(mfrow = c(2, 2))
for(sigma in seq(0.1, 0.7, len = 4)){
		y.sim <- matrix(0, Nsim, Nx)
	for(j in 1:Nsim){
		w.noise <- w + rnorm(Nx, mean = 0, sd = 0.01)
		y.sim[j, 1] <- runif(1, -5, 5)
		c <- runif(1, 0, 10) 
		y.sim[j, 2] <- y.sim[j, 1] + abs(rnorm(1, mean = 0, sd = 0.01))
		for (i in 3 : Nx){
		 	y.sim[j, i] <- y.sim[j, i-1] + (y.sim[j, i-1] - y.sim[j, i-2]) * exp(w.noise[i] * delta - sigma^2/2 * delta + sigma * rnorm(1, 0, sqrt(delta))) 
	    }
	}
	
	plot(xx, y.sim[1, ], type = "l", ylim = c(-5, 15), main = paste("sigma =", sigma) )
	for(i in 2 : Nsim){
		points(xx, y.sim[i, ], type = "l", lty = i, col = i)
	}

}
############################################################################################### 
## Fit ODE model - Basis Expansion
# list3d <- rep(0, Nt*N.auc *2)
# list3d <- array(list3d, c(Nt, N.auc), 2))
# list3d[, , 1] <- Acc.smooth
# list3d[, , 2] <- Velo.smooth
# dimnames(list3d)[[1]] <- plotpoints
# time.range <- c(0, 7)
# time.basis <- create.fourier.basis(time.range, nbasis = 21)
# lfd <- vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval = time.range)
# time.lfd<- smooth.basisPar(plotpoints,  list3d , time.basis, Lfdobj = lfd, lambda = 0.01)$fd
# Acc.fd <- time.lfd[, 1]
# Velo.fd <- time.lfd[, 2]
# fit.fd <- fRegress(Acc.fd ~ Velo.fd)

# const <- rep(0, dim(Acc.fd$coef)[2])
# xfdlist <- list(const = const, velofd = Velo.fd)
# beta0 <- with(Acc.fd, fd(basisobj = time.basis, fdname = NULL))
# beta1 <- with(Velo.fd, fd(basisobj = time.basis, fdnames = NULL))
# betalist <- list(const = fdPar(beta0), velofd = fdPar(beta1))
# fit.fd <- fRegress(Acc.fd, xfdlist, betalist)  

# ## Problem!!! This function could not eliminate intercept
# plot(plotpoints, eval.fd(plotpoints, fit.fd$betaestlis$Velo.fd$fd))
######################################################################
## Fit ODE model -- Point Wise
weight0 <- rep(0, Nt)
rsq <- rep(0, Nt)
for( i in 1: Nt){
	Y <- Acc.smooth[, i]
	X <- Velo.smooth[, i]
	weight0[i] <- as.numeric(lm(Y~X - 1)$coef)
	rsq[i] <- summary(lm(Y ~ X - 1 ))$r.square
} 
plot(plotpoints, weight0, type = "l")

## Graphicallly check fit
par(mfrow = c(1, 2))
plot(plotpoints, Acc.smooth[1, ], ylim = c(-5, 3), type = "l")
for(i in 2 : dim(Acc.smooth)[1]){
	points(plotpoints, Acc.smooth[i, ], type = "l", lty = i, col = i)
}
plot(plotpoints, Velo.smooth[1, ] * weight0, ylim = c(-5, 3), type = "l")
for(i in 2 : dim(Velo.smooth)[1]){
	points(plotpoints, Acc.smooth[i, ] - Velo.smooth[i, ] * weight0, type = "l", lty = i, col = i)
}
par(mfrow = c(1, 1))
plot(plotpoints, rsq, type = "l")

Fratio <- rep(0, Nt)
for(i in 1:Nt){
	psse0 <- sum(Acc.smooth[, i] ^ 2)
	psseL <- sum((Velo.smooth[, i] * weight0[i] - Acc.smooth[, i]) ^ 2)
	Fratio[i] <- (psse0-psseL)/psse0 
}
plot(plotpoints, Fratio, type = "l")

##################################################################################
## Estimate Stochastic Model
## d(Velo) = mu(t) * Velo dt + sigma(t) * Velo dWt

## Local kernel function, h is the max memeory considered
kern <- function(h, u){(h^2 - u^2)*(u>-h && u <=0)}

## add more evaluation points locally --> divide each interval into kk sub-intervals
kk <- 10
addpoints <- seq(0,7,len=kk*(Nt-1)+1)
LogPrice.add <-  matrix(0, N.auc, length(addpoints))
Velo.add <-  matrix(0, N.auc, length(addpoints))
Acc.add <-  matrix(0, N.auc, length(addpoints))

for(i in 1 : N.auc){
  temp <- SPLINE(EBAY[[i]]$Time,log(EBAY[[i]]$Price), addpoints)
  LogPrice.add[i, ] <- temp[, 1]
  Velo.add[i, ] <- temp[, 2]
  Acc.add[i, ] <- temp[, 3]
}

## Local weighted regression
mu <- rep(0, length(plotpoints))
rsq <- rep(0, length(plotpoints))

## max number of history retrieved: 
kk.back <- kk*20
for(i in 1:length(plotpoints)){
	response <- Acc.add[ , max(1, kk*i+1-kk.back) : (kk*i+1)]
	x <- Velo.add[, max(1, kk*i+1-kk.back) : (kk*i+1)]
	response <- as.vector(t(response))
	x <- as.vector(t(x))
	timediff <- seq(-length(x)/N.auc + 1 , 0) * (7 / ((Nt-1) * kk))
	weight <- rep(0, length(timediff))
	for(j in 1:length(weight)){
		weight[j] <- kern(7/Nt, timediff[j])
	}


	weight  <- rep( weight , N.auc )
	fit <- lm(response ~ x - 1, weights = weight)
	mu[i] <- fit$coefficient
	rsq[i] <- summary(fit)$r.squared	
}

plot(plotpoints, mu, type = "l")
plot(plotpoints, rsq, type = "l")

## Estimate volitility
## Et = (Y - Yhat) / Xt / sqrt( delta t)  --> should be normal

Et <- matrix(0, N.auc, Nt)
sigma <- rep(0, Nt)
sigma2 <- rep(0, Nt)
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

remove_quantile <- function(x, first, last){
	qnt <- quantile(x, probs=c(first, last))
	y <- x
	y[x < qnt[1]] <- NA
	y[x > qnt[2]] <- NA
	y
}

table <- matrix(0, 190, Nt)
for(i in 1:length(plotpoints)){
	Et[ , i] <- (Acc.smooth[, i] - mu[i]*Velo.smooth[, i])/Velo.smooth[,i]
	## remove real outlier?
	#temp <- remove_outliers(Et[ , i])
	## or remove some percentile
	temp <- remove_quantile(Et[, i], 0.2, 0.8)
	table[which(is.na(temp)) ,i] <- 1 
	temp2 <- Et[, i] * Velo.smooth[,i]
	sigma[i] <- sd(temp[which(!is.na(temp))] / sqrt(7/Nt))
	sigma2[i] <- sd(temp2) / sqrt(7/Nt)
}	
par(mfrow=c(1,3))
plot(plotpoints, sigma,type = "l")
plot(plotpoints, apply(table, 2,sum),type = "l",main= "no of outliers deleted")
plot(plotpoints, sigma2,type = "l")

Fratio.sto1 <- rep(0, Nt)
Fratio.sto2 <- rep(0, Nt)
Fratio.C <- rep(0, Nt)
psseL2 <- psseL1 <- 0
for(i in 1:Nt){
	psse0 <- sum(Acc.smooth[, i] ^ 2)
	psseC <- sum((Velo.smooth[, i] * mu[i]- Acc.smooth[, i]) ^ 2)
	psseL2 <- psseL1 <- 0
	psseL2 <- sum((Velo.smooth[, i] * mu[i] + sigma2[i] * rnorm(N.auc,0,sqrt(7/Nt)) - Acc.smooth[, i]) ^ 2)

	
	psseL1 <- sum((Velo.smooth[, i] * mu[i] + sigma[i] * Velo.smooth[,i] * rnorm(N.auc,0,sqrt(7/Nt)) - Acc.smooth[, i]) ^ 2)
	Fratio.C[i] <- (psse0-psseC)/psse0
	Fratio.sto1[i] <- (psse0-psseL1)/psse0 
	Fratio.sto2[i] <- (psse0-psseL2)/psse0
}
plot(plotpoints, Fratio.sto2, type = "l")
lines(plotpoints, Fratio.sto1, type = "l", col = "red")
lines(plotpoints, Fratio.C,type = "l", col = "blue", lty = 2)
lines(plotpoints, Fratio, type= "l",lty=2)



####################################################################################
## Function to compare actual values with fitted curves
## Note 7 day length is hard coded here
## Note here it should not be Acc!!! when used it should be the fitted price
compReal <- function(Acc, plotpoints, EBAY){
  if(dim(Acc)[1] != length(EBAY)){
    print("Length not match--Ebay")
    return
  }
  if(dim(Acc)[2] != length(plotpoints)){
    print("Length not match--plotpoints")
    return
  }
  # compare case by case, loop through every auction
  # save rsq to rsq[i]
  timediff <- 7/(length(plotpoints)-1)
  rsq <- rep(0, length(EBAY))

  for(i in 1:length(EBAY)){
    time <- EBAY[[i]]$Time
    logprice <- log(EBAY[[i]]$Price)
    logprice2 <- Acc[i,]
    
    for(j in 1:length(logprice)){
      # lower count of plotpoints for time[j]
      index <- as.integer( time[j]/timediff ) 
      # actual value of the plotpoints
      a <- index * 7/length(plotpoints)
      # linear estimation of fitted value
      # note index starts from 0, so need to add another 1.
      fitvalue <- (time[j]-a)/(timediff) * (logprice2[index+2] - logprice2[index+1]) + logprice2[index+1]
      rsq[i] <- rsq[i]+(fitvalue-logprice[j])^2 
    }
    
  }
  return(rsq)
  
}

## Check original fit
errorlist <- compReal(LogPrice.smooth, plotpoints, EBAY)
plot(errorlist)
sum(errorlist)

## check ODE fit
logprice.ode <- matrix(0, length(EBAY), length(plotpoints))
delta <- 7/ (length(plotpoints) - 1 )
for(i in 1:length(EBAY)){
  logprice.ode[i, 1] <- log(EBAY[[i]]$Price[1])
  logprice.ode[i, 2] <- logprice.ode[i, 1] + Velo.smooth[i,1]*delta
  for (j in 3 : length(plotpoints)){
    logprice.ode[i, j] <- (weight0[j] * delta * logprice.ode[i, j-1] - 2 * logprice.ode[i, j-1] + logprice.ode[i, j-2]) / (weight0[j] * delta -1)
  }
}
plot(plotpoints, logprice.ode[1, ], type = "l", ylim = c(-5,15))
for(i in 1:length(EBAY)){
  points(plotpoints, logprice.ode[i, ], type="l", col = i)
  i = i+1
}
error.ode <- compReal(logprice.ode, plotpoints, EBAY)
plot(error.ode)
sum(error.ode)


## check SDE fit
logprice.sde <- matrix(0, length(EBAY), length(plotpoints))
delta <- 7/ (length(plotpoints) - 1 )
#sigma0 = c(rep(0.1,200),rep(0,1400),rep(0.3,400))
rep <- 5
error.sde<-matrix(0, rep, length(EBAY))
velo.sde <- matrix(0,length(EBAY),length(plotpoints) )
for(count in 1:rep){
  for(i in 1:length(EBAY)){
    logprice.sde[i, 1] <-  log(EBAY[[i]]$Price[1])
    logprice.sde[i, 2] <- logprice.sde[i, 1] + Velo.smooth[i,1]*delta
    for (j in 3 : length(plotpoints)){
      logprice.sde[i, j] <- logprice.sde[i, j-1] + (logprice.sde[i, j-1] - logprice.sde[i, j-2]) * exp(mu[j] * delta - sigma[j]^2/2 * delta + sigma[j] * rnorm(1, 0, sqrt(delta))) 
    }
  }
  
  plot(plotpoints, logprice.sde[1, ], type = "l", ylim = c(-5,15))
  for(i in 1:length(EBAY)){
    points(plotpoints, logprice.sde[i, ], type="l", col = i)
  }
  error.sde[count, ] <- compReal(logprice.sde, plotpoints, EBAY)
}
error.sde2<-apply(error.sde,2, mean)
plot(error.sde2)
sum(error.sde2)

######## Compare with what the plots should be like!
plot(EBAY[[1]]$Time, log(EBAY[[1]]$Price), type="l", ylim = c(-5,15))
for(i in 1:190){
  points(EBAY[[i]]$Time, log(EBAY[[i]]$Price),type="l", col=i)
}