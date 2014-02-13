StepCheck <- function(time,vector,t){ 
  stage=0
  if(length(time)==1 && time[1]>=t){return(0)}
  if(length(time)==1 && time[1]<t){return(vector[1])}
  for(i in 1 : (length(time)-1)){
    if (t>=time[i] && t<time[i+1]){stage=vector[i]}
  }
  if (t>=time[length(time)]){stage=vector[length(time)]}
  return(stage)
}
#Function: Stepcheck(time, vector, t), return vector step count at t.
StepCount <- function(time,vector,t){ 
  if(length(time)==1 && time[1]>=t){return(0)}
  if(length(time)==1 && time[1]<t){return(1)}
  for(i in 1 : (length(time)-1)){
    if (t>=time[i] && t<time[i+1]){count=i}
  }
  if (t>=time[length(time)]){count=length(time)}
  return(count)
}


palm_raw <- read.csv("~/Dropbox/Ebay_2012/Modeling Online Auction/Palm+7-day+149auctions+Curve+Clustering.csv",header = TRUE)
N <- length(unique(palm_raw[,1]))
palm_id <- unique(palm_raw[, 1])
palm <- vector("list", N)
auc_length <- rep(0, N)
k <- 1
for(i in 1:dim(palm_raw)[1]){
	if(palm_raw[i, 1] == palm_id[k]){
		palm[[k]]$bid <- c(palm[[k]]$bid, palm_raw[i, 2])
		palm[[k]]$time <-c(palm[[k]]$time, palm_raw[i, 3])
	}else{
		k <- k + 1
		palm[[k]]$bid <- c(palm[[k]]$bid, palm_raw[i, 2])
		palm[[k]]$time <-c(palm[[k]]$time, palm_raw[i, 3])
		auc_length[k] <- palm_raw[i, 4] 	
	}
}

 # adding starting time/value
 for(i in 1:length(palm)){
	 palm[[i]]$time <- c(0, palm[[i]]$time)
	 palm[[i]]$bid <- c(palm[[i]]$bid[1], palm[[i]]$bid)
	
 }


 increment=function(p){
  inc=0
  if(p<=0.99){inc=0.05}
  if(p>=1 && p<=4.99){inc=0.25}
  if(p>=5 && p<=24.99){inc=0.5}
  if(p>=25 && p<=99.99){inc=1}
  if(p>=100 && p<=249.99){inc=2.5}
  if(p>=250 && p<=499.99){inc=5}
  if(p>=500 && p<=999.99){inc=10}
  if(p>=1000 && p<=2499.99){inc=25}
  if(p>=2500 && p<=4999.99){inc=50}
  if(p>=5000){inc=100}
  return(inc)
}
for( i in 1: length(palm)){
  Pricelive <- palm[[i]]$bid
  Pricelive[1] <- palm[[i]]$bid[1]
  Pricelive[2] <- palm[[i]]$bid[1]
  max <- palm[[i]]$bid[2]
  if(length(palm[[i]]$bid)>=3){
  for(j in 3: length(palm[[i]]$bid)){
    d=palm[[i]]$bid[j] - max
    if(d>0){ 
      Pricelive[j] <- min( max + increment(max), palm[[i]]$bid[j])
      max = palm[[i]]$bid[j]
      }
    if(d<0) {
      Pricelive[j] <- max(palm[[i]]$bid[j]+increment(palm[[i]]$bid[j]), Pricelive[j-1])
      }
    if(d==0){
      Pricelive[j]<-palm[[i]]$bid[j]
      }
  }
 }
 palm[[i]]$bid <- Pricelive
}
# remove bad data
#palm[[68]] <- NULL

save(palm, file = "palmAuction.rda")

N.auc = length(palm)
#################
##      See Smoothing and Simulation for some steps
#################
for(i in 1 : N.auc){
  temp <- SPLINE(palm[[i]]$time,log(palm[[i]]$bid), plotpoints)
  LogPrice.smooth[i, ] <- temp[, 1]
  Velo.smooth[i, ] <- temp[, 2]
  Acc.smooth[i, ] <- temp[, 3]
 }
par(mfrow= c(3, 4))
for(i in 1:N.auc){
	plot(plotpoints, LogPrice.smooth[i, ], col = "red", type = "l", ylim = c(-1,7))
	points(palm[[i]]$time, log(palm[[i]]$bid))
}
#########################################################
weight0 <- rep(0, Nt)
for(i in 1:Nt){
	Y <- Acc.smooth[, i]
	X <- Velo.smooth[, i]
	weight0[i] <- as.numeric(summary(lmrob(Y~X - 1))$coefficients[1,1])
}
Fratio <- rep(0, Nt)
for(i in 1:Nt){
	psse0 <- sum(Acc.smooth[, i] ^ 2)
	psseL <- sum((Velo.smooth[, i] * weight0[i] - Acc.smooth[, i]) ^ 2)
	Fratio[i] <- (psse0-psseL)/psse0 
}
plot(plotpoints, Fratio, type = "l")
