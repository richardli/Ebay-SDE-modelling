xlist <- NULL
ylist <- NULL
for(i in 1:N.auc){
	xlist[[i]] <- EBAY[[i]]$Time
	ylist[[i]] <- log(EBAY[[i]]$Price)
}


pca.spline <- function(xlist, ylist, ncomp = 3){
	knots <- c(0,1,2,3,4,5,6,6.25,6.5,6.75,6.8125,6.875,6.9375,7)
	numk<-length(knots)
	numplot <- length(plotpoints)
	norder <- 5
	nbasis <- numk + 2
	lambda <- 0.1
	numauct<-N.auc
	wbasis <- create.bspline.basis(rangeval=c(0,7),nbasis=nbasis, norder=5)
	Wfd <- list(0); length(Wfd) <- numauct
	bids <- array(0,length(knots)*length(xlist), dim = c(length(knots),length(xlist)))
	for(i in 1:length(xlist)){
		for(j in 1:length(knots)){
		  bids[j, i]=StepCheck(xlist[[i]],ylist[[i]],knots[j])	
		} 
	}
	
	
	Wfd <- Data2fd(knots,bids,basisobj=wbasis)
	coef <- as.matrix(Wfd$coefs)
	PCA <- pca.fd(Wfd, nharm = ncomp, centerfns = FALSE)
    
    # cumsum(PCA$varprop)
	# par(mfrow = c(4, 1))
	# plot.pca.fd(PCA,harm=1,pointplot=FALSE)
	# plot.pca.fd(PCA,harm=2,pointplot=FALSE)
	# plot.pca.fd(PCA,harm=3,pointplot=FALSE)
	# plot.pca.fd(PCA,harm=4,pointplot=FALSE)

	return(list(coef, PCA))
	
}


pca.all <- pca.spline(xlist, ylist, ncomp = 3)
coef.mat <- t(pca.all[[1]])
pca.ebay <- pca.all[[2]]
scores.mat <- pca.ebay$scores

cumsum(pca.all[[2]]$varprop)

feature <- cbind(scores.mat, coef.mat)

n.group <- 3
kmean.ebay <- kmeans(feature,centers = n.group)

group <- NULL
for(i in 1:n.group){
  group[[i]] <- which(kmean.ebay$cluster == i)	
}
par(mfrow = c(3,1))
for(i in 1:n.group){
	plot(xlist[[group[[i]][1]]], ylist[[group[[i]][1]]], ylim = c(-5, 5), type = "l", xlab = "time",ylab="logprice", main = paste("the ",i," group"))
	for(j in 2:length(group[[i]])){
		points(xlist[[group[[i]][j]]], ylist[[group[[i]][j]]], col = j, type ="l")
	}
	
}

subpara1 <- FitsubSDE(group[[1]])

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


FitsubSDE <- function(index){
	
	Acc.temp <- Acc.smooth[index, ]
	Velo.temp <- Velo.smooth[index, ]
    Logprice.temp <- LogPrice.smooth[index, ]
    N.auc.temp <- length(index)
    Et.temp <- matrix(0, N.auc.temp, Nt)
	sigma.temp <- rep(0, Nt)
	sigma2.temp <- rep(0, Nt)
    table <- matrix(0, N.auc.temp, Nt)
	mu.temp <- fitSDE_mu(index, kk)
	for(i in 1:length(plotpoints)){
		Et.temp[ , i] <- (Acc.temp[, i] - mu.temp[i]*Velo.temp[, i])/Velo.temp[,i]
		## remove real outlier?
		#temp <- remove_outliers(Et[ , i])
		## or remove some percentile
		temp <- remove_quantile(Et.temp[, i], 0.2, 0.8)
		table[which(is.na(temp)) ,i] <- 1 
		temp2 <- Et.temp[, i] * Velo.temp[,i]
		sigma.temp[i] <- sd(temp[which(!is.na(temp))] / sqrt(7/Nt))
		sigma2.temp[i] <- sd(temp2) / sqrt(7/Nt)
	}	
	par(mfrow=c(1,3))
	plot(plotpoints, sigma.temp,type = "l")
	plot(plotpoints, apply(table, 2,sum),type = "l",main= "no of outliers deleted")
	plot(plotpoints, sigma2.temp,type = "l")
	Fratio.sto1 <- rep(0, Nt)
	Fratio.sto2 <- rep(0, Nt)
	Fratio.C <- rep(0, Nt)
	psseL2 <- psseL1 <- 0
	for(i in 1:Nt){
		psse0 <- sum(Acc.temp[, i] ^ 2)
		psseC <- sum((Velo.temp[, i] * mu.temp[i]- Acc.temp[, i]) ^ 2)
		psseL2 <- psseL1 <- 0
		psseL2 <- sum((Velo.temp[, i] * mu.temp[i] + sigma2.temp[i] * rnorm(N.auc.temp,0,sqrt(7/Nt)) - Acc.temp[, i]) ^ 2)
		psseL1 <- sum((Velo.temp[, i] * mu.temp[i] + sigma.temp[i] * Velo.temp[,i] * rnorm(N.auc.temp,0,sqrt(7/Nt)) - Acc.temp[, i]) ^ 2)
		Fratio.C[i] <- (psse0-psseC)/psse0
		Fratio.sto1[i] <- (psse0-psseL1)/psse0 
		Fratio.sto2[i] <- (psse0-psseL2)/psse0
	}
	plot(plotpoints, Fratio.sto2, type = "l")
	lines(plotpoints, Fratio.sto1, type = "l", col = "red")
	lines(plotpoints, Fratio.C,type = "l", col = "blue", lty = 2)

return(list(sigma1 = sigma.temp, sigma2 = sigma2.temp, F0 = Fratio.C, F1 = Fratio.sto1, F2 = Fratio.sto2))
} 


fitSDE_mu <- function(index, kk){
	kern <- function(h, u){(h^2 - u^2)*(u>-h && u <=0)}
	
	addpoints <- seq(0,7,len=kk*(Nt-1)+1)
	
	Acc.temp <- Acc.smooth[index, ]
	Velo.temp <- Velo.smooth[index, ]
    Logprice.temp <- LogPrice.smooth[index, ]
    N.auc.temp <- length(index)
	LogPrice.add <-  matrix(0, N.auc.temp, length(addpoints))
	Velo.add <-  matrix(0, N.auc.temp, length(addpoints))
	Acc.add <-  matrix(0, N.auc.temp, length(addpoints))
	
	for(i in 1 : N.auc.temp){
	  temp <- SPLINE(EBAY[[index[i]]]$Time,log(EBAY[[index[i]]]$Price), addpoints)
	  LogPrice.add[i, ] <- temp[, 1]
	  Velo.add[i, ] <- temp[, 2]
	  Acc.add[i, ] <- temp[, 3]
	}
	
	mu <- rep(0, length(plotpoints))
    rsq <- rep(0, length(plotpoints))

	kk.back <- kk*10
	for(i in 1:length(plotpoints)){
		response <- Acc.add[ , max(1, kk*(i-1)+1-kk.back) : (kk*(i-1)+1)]
		x <- Velo.add[, max(1, kk*(i-1)+1-kk.back) : (kk*(i-1)+1)]
		response <- as.vector(t(response))
		x <- as.vector(t(x))
		timediff <- seq(-length(x)/N.auc.temp + 1 , 0) * (7 / ((Nt-1) * kk))
		weight <- rep(0, length(timediff))
		for(j in 1:length(weight)){
			weight[j] <- kern(7/Nt, timediff[j])
		}
	
	
		weight  <- rep( weight , N.auc.temp )
		fit <- lm(response ~ x - 1, weights = weight)
		mu[i] <- fit$coefficient
		rsq[i] <- summary(fit)$r.squared	
	}
	
	return(mu)

}


