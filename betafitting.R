############################################## Include beta fitting and regression##########################
#################################### Parameters: START, NT.grid, EBAY.training##############################
install.packages(c("optimx","BB","ucminf","Rcgmin","quadprog","Rvmmin","minqa"))
library(optimx)
library(BB)
library(ucminf)
library(Rcgmin)
library(quadprog)
library(Rvmmin)
library(minqa)
######################
START=0
K=2
NT=500
END=7  #For the time being, the end could not be changed!!!
## Beta fitting
######################
smoothF.1=NULL
pdf.1.beta=NULL
acc <- NULL
Fpar=matrix(0,length(palm_sev),2)
par(mfrow = c(3,4))
for(i in 1:length(palm_sev)){
  
{  ## Obtain the parameter
  
  lastday=StepCount(palm_sev[[i]]$time,palm_sev[[i]]$time, START)
  if(lastday==length(palm_sev[[i]]$time)){
    Fpar[i,]=c(9999,9999)
    smoothF.1[[i]]=rep(0,NT)
    pdf.1.beta[[i]]=rep(0,NT)
    next
  }
  
  #  To handle the last term being inf.
  deltaT=max(0.000001,(palm_sev[[i]]$time[length(palm_sev[[i]]$time)]-palm_sev[[i]]$time[length(palm_sev[[i]]$time)-1]))
  x1= min(palm_sev[[i]]$time[length(palm_sev[[i]]$time)]+deltaT,END)
  xx= c(palm_sev[[i]]$time[lastday:length(palm_sev[[i]]$time)],x1)
  yy= c(log(palm_sev[[i]]$bid)[lastday:length(log(palm_sev[[i]]$bid))])
  y1= 2*yy[length(yy)]-yy[length(yy)-1]
  yy=c(yy,y1)
  x= (xx-START)/(END-START)
  y= (yy-min(yy))/(max(yy)-min(yy))
  x[1]=0
  
  
  DistB=function(a){
    dist=sum((y-pbeta(x,a[1],a[2]))^2)+sum((x-qbeta(y,a[1],a[2]))^2)
  }  
  DistA=function(a){
    (mean(x)-a[1]/(a[1]+a[2]))^2+(var(x)-a[1]*a[2]/((a[1]+a[2])^2*(a[1]+a[2]+1)))^2 
  }
  
  para0 =unlist(optimx(c(1,0), DistA, method="BFGS")$par)
  para0=abs(para0)
  para = optimx(para0, DistB , method="BFGS")
  if (para$conv == 9999){para = optimx(para0, DistB , method="Nelder-Mead")}
  Fpar[i,] = unlist(para$par)
}


smoothF.1.temp=NULL
pdf.1.beta.temp=NULL
acc.smooth.temp <- NULL
for(j in 1:NT){
  integrand=function(u){u^(Fpar[i,1]-1)*(1-u)^(Fpar[i,2]-1)}
  smoothF.1.temp=cbind(smoothF.1.temp,integrate(integrand, 0, j/(NT+1))$value/beta(Fpar[i,1],Fpar[i,2]))
  pdf.1.beta.temp = cbind(pdf.1.beta.temp, integrand(j/(NT+1))/beta(Fpar[i,1],Fpar[i,2]))
	acc.smooth.temp <- cbind(acc.smooth.temp, ((Fpar[i,1] - 1) * (j/(NT+1))^(Fpar[i, 1] -2)*(1-(j/(NT+1)))^(Fpar[i,2]-1) - (j/(NT+1))^(Fpar[i,1] - 1) * (Fpar[i,2]-1)*(1-(j/(NT+1)))^(Fpar[i, 2]-2))/beta(Fpar[i,1],Fpar[i,2]))


}





# smoothF.1[[i]]=smoothF.1.temp * (max(yy)-min(yy)) + min(yy)
# pdf.1.beta[[i]]=pdf.1.beta.temp * (max(yy)-min(yy))/(END-START)

## Change to last day, edited Jan 2013
smoothF.1[[i]]<-(smoothF.1.temp * (max(yy)-min(yy)) + min(yy))#[4001:8000]
pdf.1.beta[[i]]<-(pdf.1.beta.temp * (max(yy)-min(yy))/(END-START))#[4001:8000]
acc[[i]] <- (acc.smooth.temp *(max(yy)-min(yy))/(END-START)) 
#par(mfrow=c(1,2))
xaxis=seq(START,END,len=NT)
plot(xaxis, smoothF.1[[i]],type="l",xlim=c(START,END))
points(palm_sev[[i]]$time,log(palm_sev[[i]]$bid),cex=0.5)
points(x*(END-START)+START,y* (max(yy)-min(yy)) + min(yy),col="red")
plot(xaxis, pdf.1.beta[[i]],type = "l")

#Print out the plots of Fitting and lambda
#  par(mfrow=c(1,2))
#  xaxis=seq(START,END,len=NT)
#  plot(xaxis, smoothF.1[[i]],type="l",xlim=c(START,END),ylim=c(0,1))
#  points(palm_sev[[i]]$time,palm_sev[[i]]$bid,cex=0.5)
#  points(x*(END-START)+START,y* (max(yy)-min(yy)) + min(yy),col="red")
#  plot(xaxis, pdf.1.beta[[i]],cex=0.4)

}

N.auc <- length(palm_sev)
logprice.smooth <- matrix(0, N.auc, NT)
velo.smooth <- matrix(0, N.auc, NT)
acc.smooth <- matrix(0, N.auc, NT)


for(i in 1:N.auc){
	logprice.smooth[i, ] <- smoothF.1[[i]]
	velo.smooth[i, ] <- pdf.1.beta[[i]]
	acc.smooth[i, ] <- acc[[i]]
}
