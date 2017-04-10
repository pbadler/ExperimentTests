rm(list=ls(all=TRUE)); graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

# read in distance weights
dists <- read.csv(paste0(root,"/ExperimentTests/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))

doSpp <- "PSSP"
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste(root,"/ExperimentTests/data/idaho",sep="")
dataDir2 <- paste(root,"/ExperimentTests/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

# import old data--------------------------------------------------------

setwd("growth")
source("fetchGrowthData.r")

D1 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"

############## Fit a simple growth model. 
g1 <- lm(logarea.t1 ~ logarea.t0 + factor(Group) + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts,data=D1); 

############## Generate data that satisfy the fitted growth model (ignoring Group and competition) 
x = rnorm(nrow(D1), mean = mean(D1$logarea.t0),sd=sd(D1$logarea.t0)); 
yhat = g1$coef[1]+g1$coef[2]*x; y = yhat + rnorm(length(x),mean=0,sd=sqrt(mean(g1$residuals^2))); 

############## Two ways of looking for bias in the fitted model 
dev.new(height=11,width=8); 
par(mfcol=c(2,1),mar=c(4,4,1,1),mgp=c(2,1,0)); 

# (1) Plot predicted (y) versus observed (x) 
scatter.smooth(y,yhat,xlab="Observed log size", ylab="Predicted log size"); 
abline(0,1,lty=2,col="red",lwd=2); 
legend("topleft",legend=c("Trend line", "1:1 line"),col=c("black","red"),lty=c(1,2),lwd=c(1,2)); 

# (2) Plot observed (y) versus predicted (x) 
scatter.smooth(yhat,y,xlab="Predicted log size", ylab="Observed log size"); 
abline(0,1,lty=2,col="red",lwd=2); 
legend("topleft",legend=c("Trend line", "1:1 line"),col=c("black","red"),lty=c(1,2),lwd=c(1,2)); 

title(main="PSSP historical data"); 
dev.copy2pdf(file="GrowthBiasPlot.pdf"); 

##############################################
# Effect of the discrete size classes 
##############################################

graphics.off(); 
dev.new(width=9,height=9); par(mfrow=c(2,2),mgp=c(2,1,0),cex.axis=1.3,cex.lab=1.3); 

for(ispp in 1:4) {

doSpp <- sppList[ispp]; 
D1 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)

# Fit a simple linear growth model 
g0 <- lm(logarea.t1 ~ logarea.t0 + factor(Group) + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts,data=D1); 

# Quadratic growth model model 
g1 <- lm(logarea.t1 ~ logarea.t0 + I(logarea.t0^2) + factor(Group) + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts,data=D1); 

# Quadratic growth model omitting discrete small size classes 
e=which(D1$logarea.t0 > 0.5); 
D2 = D1[e,]; 
g2 <- lm(logarea.t1 ~ logarea.t0 + I(logarea.t0^2) + factor(Group) + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts,data=D2); 

plot(logarea.t1 ~ logarea.t0,data=D1,xlab="Log area(t)", ylab="Log area(t+1)");  

x = sort(D1$logarea.t0); 
y1 = g1$coef[1] + g1$coef[2]*x + g1$coef[3]*x^2; 
y2 = g2$coef[1] + g2$coef[2]*x + g2$coef[3]*x^2; 
matpoints(x,cbind(y1,y2),type="l",col=c("black","red"),lwd=2); 

abline(g0$coef[1],g0$coef[2],col="blue",lty=2,lwd=2); 

if(ispp==1) legend("topleft",legend=c("Linear, all", "Quadratic, all", "Quadratic, Area(t)>0.5"),col=c("blue","black","red"),lwd=2,bty="n"); 

# Estimate Jensen factor for large plants 
A = quantile(D2$logarea.t0,0.5); 
e=which(D1$logarea.t0 > A); 
D3 = D1[e,]; 
g3 <- lm(logarea.t1 ~ logarea.t0 + I(logarea.t0^2) + factor(Group) + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts,data=D3); 

sigma2=mean(g3$residuals^2); JenFac=exp(sigma2/2); 
title(main=paste(doSpp,"historical data, J=",round(JenFac,digits=2))); 

}

dev.copy2pdf(file="GrowthModelCheck.pdf"); 