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
	PCA <- pca.fd(Wfd, nharm = ncomp)
    
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

feature <- cbind(scores.mat, coef.mat)

n.group <- 3
kmean.ebay <- kmeans(feature,centers = n.group)

group <- NULL
for(i in 1:n.group){
  group[[i]] <- which(kmean.ebay$cluster == i)	
}
par(mfrow = c(2,1))
for(i in 1:n.group){
	plot(xlist[[group[[i]][1]]], ylist[[group[[i]][1]]], ylim = c(-5, 5), type = "l", xlab = "time",ylab="logprice", main = paste("the ",i," group"))
	for(j in 2:length(group[[i]])){
		points(xlist[[group[[i]][j]]], ylist[[group[[i]][j]]], col = j, type ="l")
	}
	
}





