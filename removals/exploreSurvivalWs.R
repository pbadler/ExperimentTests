# PBA March 2016

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

#####################
#  1. Calculate W's
####################

sppList <- c("ARTR","HECO","POSE","PSSP")
dataDir <- paste(root,"/driversdata/data/idaho",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

allD <- list(NULL)

for(j in 1:length(sppList)){

  # import survival data
  doSpp <- sppList[j]
  survDfile=paste(dataDir,"/speciesdata/",doSpp,"/survD.csv",sep="")
  survD=read.csv(file=survDfile)
  D=survD[survD$allEdge==0,];
  D$year <- D$year
  D$logarea=log(D$area)
  D$quad=as.character(D$quad)
  
  # import neighbor data
  ringD <- read.csv(paste(dataDir,"/speciesdata/",doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
  ringD$year<-ringD$year
  
  # merge D with ringD (D contains fewer rows)
  D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
  D=D[order(D$X),]
  rm(ringD,survD)
  row.names(D) <- NULL  
  
  # calculate W's (MAKE SURE NOT TO REORDER D!)
  W <- matrix(NA,NROW(D),length(sppList))
  colnames(W) <- paste("W.",sppList,sep="")
  dists <- read.csv(paste(dataDir,"/speciesdata/IdahoDistanceWeights.csv",sep=""));
  for(iSpp in 1:length(sppList)){
    neighborCols=which(substr(names(D),1,4)==sppList[iSpp]) # pull out annulus data for the focal species 
    dist_wts<- dists[,paste0(sppList[iSpp])]
    C <- data.matrix(D[,neighborCols]) #matrix of conspecific areas in the annuli 
    W[,iSpp] <- C%*%dist_wts 
  }
  
  # reformat D
  D <- D[,c("X","quad","year","trackID","area","survives","age","distEdgeMin","allEdge","seedling","QuadName","Grazing","Group","logarea","species")]
  D <- cbind(D,W)
  
  allD[[j]] <- D
  
} # next j species

names(allD) <- sppList

#####################
#  2. Visualize W's
####################

Wcols <- grep("W.*",names(allD[[1]]))

# histograms
pdf("W-histograms.pdf",height=6,width=6)
par(mfcol=c(2,4),mgp=c(2,0.5,0),tcl=-0.2,oma=c(0,0,2,0))
for(i in 1:length(sppList)){
  for(j in 1:length(sppList)){
    hist(allD[[i]][,Wcols[j]],breaks=20,xlab="W",main=names(allD[[i]])[Wcols[j]])
    hist(sqrt(allD[[i]][,Wcols[j]]),breaks=20,xlab="sqrt(W)",main="")
  }
  mtext(sppList[i],side=3,outer=T)
}
dev.off()

# scatter plot matrix
pdf("W-scatter.pdf",height=6,width=6)
for(i in 1:length(sppList)){
  pairs(sqrt(allD[[i]][,Wcols]),main=sppList[i])
}
dev.off()

# wind rose (four dimensions)
# first format data for segments()
Nspp <-length(sppList)
xyDat <- matrix(0,NROW(D),2*Nspp + 2)
xyDat[,2] <- D[,Wcols[1]]
xyDat[,3] <- D[,Wcols[2]]
xyDat[,6] <- D[,Wcols[3]]
xyDat[,7] <- D[,Wcols[4]]
xyDat[,10] <- D[,Wcols[1]]
xyDat <- sqrt(xyDat)
xyDat[,6] <- -1*xyDat[,6]; xyDat[,7] <- -1*xyDat[,7]
xyLong <- matrix(as.vector(t(xyDat)),nrow=NROW(xyDat)*5,2,byrow=T)

# plot on sqrt scale
maxW <- max(abs(xyLong))*1.02
myTics <- c(-round(maxW),-round(maxW/2),round(maxW),round(maxW/2))
par(mfrow=c(1,1),mgp=c(2,0.5,0),tcl=-0.2)
plot(x=0,y=0,xlim=c(-1*maxW,maxW),ylim=c(-1*maxW,maxW),type="n",main=doSpp,axes=F,xlab="",ylab="")
axis(1,pos=0,at=myTics); axis(2,pos=0,at=myTics);
mtext("sqrt(W.ARTR)",side=3,at=0.5)
mtext("sqrt(W.HECO)",side=4,at=0.5)
mtext("sqrt(W.POSE)",side=1,at=0.5)
mtext("sqrt(W.PSSP)",side=2,at=0.5)
for(i in 1:NROW(D)){
  lines(xyLong[(1+(i-1)*5):(5+(i-1)*5),],col="#0000FF07",lwd=2)
}



