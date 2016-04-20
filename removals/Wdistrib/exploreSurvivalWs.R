# PBA March 2016

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/Wdistrib",sep="")); # modify as needed 

#########################################
#  1. Import data and calculate W's
#########################################

sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste(root,"/driversdata/data/idaho",sep="")
dataDir2 <- paste(root,"/driversdata/data/idaho_modern",sep="")

# set up distance weights
dists <- read.csv(paste(dataDir1,"/speciesdata/IdahoDistanceWeights.csv",sep=""));
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

# import old data
setwd("..")
source("survival/fetchSurvData.r")
setwd("Wdistrib")

allD <-list()
for(iSpp in 1:4){
  
  doSpp <- sppList[iSpp]
  
  # get old data
  D1 <- fetchSdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
  D1$year <- D1$year+1900
  
  # import modern data
  D2 <- fetchSdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)
  
  # combine old and modern
  Dall <- rbind(D1,D2)
  
  # merge in treatment data
  tmp <- read.csv(paste(dataDir2,"/quad_info.csv",sep=""))
  tmp <- tmp[,c("quad","Treatment","Group")]
  Dall <- merge(Dall,tmp, all.x=T)
  
  # account for removal in baseline years
  #   ii <- which(Dall$year>=2011 & Dall$Treatment=="No_shrub")
  #   Dall$W.ARTR[ii] <- 0
  #   ii <- which(Dall$year>=2011 & Dall$Treatment=="No_grass")
  #   Dall$W.HECO[ii] <- 0 ; Dall$W.POSE[ii] <- 0 ; Dall$W.PSSP[ii] <- 0

  allD[[iSpp]] <- Dall[,c("year","Treatment","Group","W.ARTR", "W.HECO","W.POSE","W.PSSP","W.allcov","W.allpts")]
  rm(Dall,D1,D2,tmp)

}

names(allD) <- sppList[1:4]

#####################
#  2. Model W's (with help from Giles Hooker)
####################

doSpp <- 4 # look at PSSP

# grab control treatment data
cD <- subset(allD[[doSpp]],Treatment=="Control") 
groupI <- which(cD$Group=="E1") # keep track of this for later
Wdat = cD[,4:ncol(cD)]
Wnonzero = Wdat > 0

#scatter plots
pairs(Wdat)

# grab data from no shrub plots
trtD <- subset(allD[[doSpp]],Treatment=="No_shrub" & year==2011)  
Wdat.trt = trtD[,4:ncol(trtD)]
Wnonzero.trt = Wdat.trt > 0
pairs(Wdat.trt)

myCol <- c(rep("black",dim(Wdat)[1]),rep("red",dim(Wdat.trt)[1]))
allWdat <-rbind(Wdat,Wdat.trt)
X11()
pairs(allWdat,col=myCol)



# Brief analysis -- are there correlations among whether there is a zero or not
cor(Wnonzero) # All pretty tiny. 
cor(Wnonzero.trt) # pretty strong negative correlation with HECO

# To analyze the rest, we'll look at Box-Cox transformations - these find lambda 
# so that X^lambda is as close to normal as possible. boxcox() in R just produces
# a list of lambdas and the corresponding likelihood so we have to do a bit of
# processing
library(MASS)
W.bc = Wdat  # will store transformed data
W.bc.trt = Wdat.trt
for(i in 1:6){
  
  # Select non-zero entries for this column
  t.dat = Wdat[Wnonzero[,i],i]
  
  # boxcox
  bc = boxcox(t.dat~1,lambda = seq(-2,2,by=0.01),plotit=FALSE)
 
  # optimal lambda
  lambda = bc$x[which.max(bc$y)]
  
  # Add into W.bc
  W.bc[Wnonzero[,i],i] = t.dat^lambda
  
  # apply same transformation to treatment plots
  W.bc.trt[Wnonzero.trt[,i],i] = Wdat.trt[Wnonzero.trt[,i],i]^lambda
  
}

# compare distributions of control and treatment plot W's

filename <- paste0(sppList[doSpp],"_W_byTrt.pdf")
pdf(filename,height=5,width=8.5)
  par(mfrow=c(2,3),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0))
  for(i in 1:6){
    
    pNonZero <-sum(Wnonzero[groupI,i])/length(Wnonzero[groupI,i]) # just E1 control quads
    #pNonZero <-sum(Wnonzero[,i])/length(Wnonzero[,i])  # all control quads
    pNonZero.trt <- sum(Wnonzero.trt[,i])/length(Wnonzero.trt[,i]) 
    dens <- density(W.bc[which(Wnonzero[,i]==1 & cD$Group=="E1"),i]) # just E1 control quads
    #dens <- density(W.bc[Wnonzero[,i],i])  # all control quads
    dens.trt <- density(W.bc.trt[Wnonzero.trt[,i],i])
    ylims <- c(0,1.2*max(c(dens$y,dens.trt$y)))
    xlims <- c(min(c(dens$x,dens.trt$x)),max(c(dens$x,dens.trt$x)))
    plot(dens,xlab="",ylab="",,xlim=xlims,ylim=ylims,main=names(cD)[3+i],col="black")
    lines(dens.trt,col="red")
    legend("topleft",c(paste0("Control, p(W>0)=",round(pNonZero,2)),
                       paste0("No shrub, p(W>0)=",round(pNonZero.trt,2))),
           lty=1,col=c("black","red"),bty="n")
    
  }
  
  mtext("Box-Cox W",side=1,line=0.5,outer=T)
  mtext("Density",side=2,line=0.5,outer=T)

dev.off()

# A first analysis: correlation based on rows where there are no zeros
W.all = apply(Wnonzero,1,prod)
cor(W.bc[W.all==1,])
# Just a few correlations larger than 0.1, none over 0.15

# do it again with just PSSP and ARTR
W.all = apply(Wnonzero[,c(1,4)],1,prod)
cor(W.bc[W.all==1,c(1,4)])
# almost 0.19

# Further analyses are based on conditionals; here we'll regress each column on 
# the others, allowing a separate effect for being zero in each. We'll also 
# predict whether or not the response is non-zero via a logistic regression. 

lm.mods = list()       # Models for non-zeros
glm.mods = list()      # Models to predict whether or not you're zero
Rsq = rep(0,4)         # R-squared for non-zero model
ave.dev = rep(0,4)     # Average drop in deviance per observation

for(i in 1:6){
 # Linear regression for non-zeros
 t.dat = data.frame(y = W.bc[Wnonzero[,i],i],W.bc[Wnonzero[,i],-i], Wnonzero[Wnonzero[,i],-i])
 
 lm.mods[[i]] = lm(y~.,data=t.dat)

 print(summary(lm.mods[[i]]))
 Rsq[i] = summary(lm.mods[[i]])$r.squared

 # Generalized Linear regression for whether or not the record is zero
 t.dat = data.frame(y = Wnonzero[,i],W.bc[,-i], Wnonzero[,-i])
 glm.mods[[i]] = glm(y~.,data=t.dat,family='binomial')

 print(summary(glm.mods[[i]])) 
 
 ave.dev[i] = (glm.mods[[i]]$null.deviance - glm.mods[[i]]$deviance)/nrow(t.dat)
}

# Some of these effects are significant, but the change in explanatory power is 
# tiny


#####################
#  3. Visualize W's
####################

Wcols <- grep("W.*",names(allD[[1]]))

# histograms
pdf("W-histograms.pdf",height=10,width=4)
par(mfrow=c(6,2),mgp=c(2,0.5,0),tcl=-0.2,oma=c(0,0,2,0))
for(i in 1:4){
  for(j in 1:length(sppList)){
    hist(allD[[i]][,Wcols[j]],breaks=20,xlab="W",main=names(allD[[i]])[Wcols[j]])
    hist(sqrt(allD[[i]][,Wcols[j]]),breaks=20,xlab="sqrt(W)",main="")
  }
  mtext(sppList[i],side=3,outer=T)
}
dev.off()

# scatter plot matrix
pdf("W-scatter.pdf",height=6,width=6)
for(i in 1:4){
  pairs(sqrt(allD[[i]][,Wcols]),main=sppList[i])
}
dev.off()

# wind rose (four dimensions)
# myRoot <- 3
# Nspp <-length(sppList)
# pdf("W-windrose.pdf",height=2.5,width=8.5)
# par(mfrow=c(1,4),tcl=-0.2,mgp=c(2,0.5,0),mar=c(2,2,4,2))
# 
# for(i in 1:length(sppList)){
#   
#   # first format data for lines()
#   xyDat <- matrix(0,NROW(allD[[i]]),2*Nspp + 2)
#   xyDat[,2] <- allD[[i]][,Wcols[1]]
#   xyDat[,3] <- allD[[i]][,Wcols[2]]
#   xyDat[,6] <- allD[[i]][,Wcols[3]]
#   xyDat[,7] <- allD[[i]][,Wcols[4]]
#   xyDat[,10] <- allD[[i]][,Wcols[1]]
#   xyDat <- xyDat^(1/myRoot)
#   xyDat[,6] <- -1*xyDat[,6]; xyDat[,7] <- -1*xyDat[,7]
#   xyLong <- matrix(as.vector(t(xyDat)),nrow=NROW(xyDat)*5,2,byrow=T)
#   
#   # plot on sqrt scale
#   maxW <- max(abs(xyLong))*1.02
#   myTics <- c(-round(maxW),-round(maxW/2),round(maxW),round(maxW/2))
#   plot(x=0,y=0,xlim=c(-1*maxW,maxW),ylim=c(-1*maxW,maxW),type="n",main=sppList[i],axes=F,xlab="",ylab="")
#   axis(1,pos=0,at=myTics); axis(2,pos=0,at=myTics);
#   mtext("W.ARTR",side=3,at=0.5,cex=0.7)
#   mtext("W.HECO",side=4,at=0.5,cex=0.7)
#   mtext("W.POSE",side=1,at=0.5,cex=0.7)
#   mtext("W.PSSP",side=2,at=0.5,cex=0.7)
#   for(k in 1:NROW(allD[[i]])){
#     lines(xyLong[(1+(k-1)*5):(5+(k-1)*5),],col="#0000FF07",lwd=1.5)
#   }
# 
# } # next i spp
# 
# dev.off()
# 
