# PBA March 2016

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

#####################
#  1. Calculate W's
####################

# inputs
doSpp <- "POSE"
sppList <- c("ARTR","HECO","POSE","PSSP")
dataDir <- paste(root,"/driversdata/data/idaho",sep="")
survDfile=paste(dataDir,"/speciesdata/",doSpp,"/survD.csv",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# import survival data
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

#####################
#  2. Visualize W's
####################

Wcols <- grep("W.*",names(D))

# histograms
pdf("W-histograms.pdf",height=6,width=6)
par(mfcol=c(2,4),mgp=c(2,0.5,0),tcl=-0.2)
for(i in 1:length(sppList)){
  hist(D[,Wcols[i]],breaks=20,xlab="W",main=names(D)[Wcols[i]])
  hist(sqrt(D[,Wcols[i]]),breaks=20,xlab="sqrt(W)",main="")
}
dev.off()

# scatter plot matrix
pdf("W-scatter.pdf",height=6,width=6)
pairs(sqrt(D[,Wcols]),main=doSpp)
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
par(mfrow=c(1,1),mgp=c(2,0.5,0),tcl=-0.2)
plot(1,1,xlim=c(-1*maxW,maxW),ylim=c(-1*maxW,maxW),type="n",main=doSpp)
for(i in 1:NROW(D)){
  lines(xyLong[(1+(i-1)*5):(5+(i-1)*5),],col="#FF000010")
}



