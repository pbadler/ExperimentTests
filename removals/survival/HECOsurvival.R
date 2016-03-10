# PBA March 2016

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/survival",sep="")); # modify as needed 

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- "HECO"
sppList <- c("ARTR","HECO","POSE","PSSP")
dataDir1 <- paste(root,"/driversdata/data/idaho",sep="")
dataDir2 <- paste(root,"/driversdata/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# import old data--------------------------------------------------------
survDfile=paste(dataDir1,"/speciesdata/",doSpp,"/survD.csv",sep="")
survD=read.csv(file=survDfile)
D1=survD[survD$allEdge==0,];
D1$year <- D1$year
D1$logarea=log(D1$area)
D1$quad=as.character(D1$quad)

# import neighbor data
ringD <- read.csv(paste(dataDir1,"/speciesdata/",doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
ringD$year<-ringD$year

# merge D with ringD (D contains fewer rows)
D1<-merge(D1,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
D1=D1[order(D1$X),]
rm(ringD,survD)
row.names(D1) <- NULL  

# calculate W's (MAKE SURE NOT TO REORDER D!)
W <- matrix(NA,NROW(D1),length(sppList))
colnames(W) <- paste("W.",sppList,sep="")
dists <- read.csv(paste(dataDir1,"/speciesdata/IdahoDistanceWeights.csv",sep=""));
for(iSpp in 1:length(sppList)){
  neighborCols=which(substr(names(D1),1,4)==sppList[iSpp]) # pull out annulus data for the focal species 
  dist_wts<- dists[,paste0(sppList[iSpp])]
  C <- data.matrix(D1[,neighborCols]) #matrix of conspecific areas in the annuli 
  W[,iSpp] <- C%*%dist_wts 
}

# reformat D
D1 <- D1[,c("X","quad","year","trackID","area","survives","age","distEdgeMin","allEdge","seedling","QuadName","Grazing","Group","logarea","species")]
D1 <- cbind(D1,W)
D1$Treatment <- "Control"

# import modern data--------------------------------------------------------
survDfile=paste(dataDir2,"/speciesdata/",doSpp,"/survD.csv",sep="")
survD=read.csv(file=survDfile)
D2=survD[survD$allEdge==0,];
D2$year <- D2$year
D2$logarea=log(D2$area)
D2$quad=as.character(D2$quad)

# import neighbor data
ringD <- read.csv(paste(dataDir2,"/speciesdata/",doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
ringD$year<-ringD$year

# merge D with ringD (D contains fewer rows)
D2<-merge(D2,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
D2=D2[order(D2$X),]
rm(ringD,survD)
row.names(D2) <- NULL  

# calculate W's (MAKE SURE NOT TO REORDER D!)
W <- matrix(NA,NROW(D2),length(sppList))
colnames(W) <- paste("W.",sppList,sep="")
#dists <- read.csv(paste(dataDir1,"/speciesdata/IdahoDistanceWeights.csv",sep=""));
for(iSpp in 1:length(sppList)){
  neighborCols=which(substr(names(D2),1,4)==sppList[iSpp]) # pull out annulus data for the focal species 
  dist_wts<- dists[,paste0(sppList[iSpp])]
  C <- data.matrix(D2[,neighborCols]) #matrix of conspecific areas in the annuli 
  W[,iSpp] <- C%*%dist_wts 
}

# reformat D
D2 <- D2[,c("X","quad","year","trackID","area","survives","age","distEdgeMin","allEdge","seedling","QuadName","Grazing","Group","logarea","species")]
D2 <- cbind(D2,W)

# merge in treatment data
tmp <- read.csv(paste(dataDir2,"/quad_info.csv",sep=""))
tmp <- tmp[,c("quad","Treatment")]
D2 <- merge(D2,tmp, all.x=T)

# account for removal in baseline years
if(doSpp!="ARTR"){
  ii <- which(D2$year>=2011 & D2$Treatment=="No_shrub")
  D2$W.ARTR[ii] <- 0
}else{
  ii <- which(D2$year>=2011 & D2$Treatment=="No_grass")
   D2$W.HECO[ii] <- 0 ; D2$W.POSE[ii] <- 0 ; D2$W.PSSP[ii] <- 0
}

# combine old and modern
allD <- rbind(D1,D2)
rm(D1,D2,tmp,W)

# clean up dataset ----------------------------------------------
allD$year[allD$year<2000] <- allD$year[allD$year<2000] + 1900

if(doSpp=="ARTR"){
  keep <- which(is.element(allD$Treatment,c("Control","No_grass")))
}else{
  keep <- which(is.element(allD$Treatment,c("Control","No_shrub")))
}
allD <- allD[keep,]

# remove outliers?

#########################################
#  2. Fit models
#########################################

# set up indicator variables
allD$Treatment2 <- allD$Treatment
allD$Treatment2[allD$year>2000] <- "Modern"
allD$Treatment3 <- allD$Treatment
allD$Treatment3[allD$Treatment=="Control" & allD$year>2000] <- "ControlModern"

allD$year <- as.factor(allD$year)

# Fit models with INLA

library("INLA")

# Set up ID variables for INLA random effects
allD$GroupID <- as.numeric(allD$Group)
allD$yearID <- 100+as.numeric(allD$year) # for random year offset on intercept

# baseline model
m0 <- inla(survives ~ logarea+ W.ARTR + W.HECO + W.POSE + W.PSSP+
  logarea:W.ARTR +logarea:W.HECO + logarea:W.POSE + logarea:W.PSSP+
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("binomial"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

# add intercepts
m1a <- inla(survives ~ logarea+ Treatment + W.ARTR + W.HECO + W.POSE + W.PSSP+
  logarea:W.ARTR +logarea:W.HECO + logarea:W.POSE + logarea:W.PSSP+
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("binomial"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

m1b <- inla(survives ~ logarea+ Treatment2 + W.ARTR + W.HECO + W.POSE + W.PSSP+
  logarea:W.ARTR +logarea:W.HECO + logarea:W.POSE + logarea:W.PSSP+
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("binomial"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

m1c <- inla(survives ~ logarea+ Treatment3 + W.ARTR + W.HECO + W.POSE + W.PSSP+
  logarea:W.ARTR +logarea:W.HECO + logarea:W.POSE + logarea:W.PSSP+
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("binomial"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

# add treatment intercepts and interactions with size
m2a <- inla(survives ~ logarea+ Treatment + logarea:Treatment+ W.ARTR + W.HECO + W.POSE + W.PSSP+
  logarea:W.ARTR +logarea:W.HECO + logarea:W.POSE + logarea:W.PSSP+
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("binomial"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

m2b <- inla(survives ~ logarea+ Treatment2 + logarea:Treatment2+ W.ARTR + W.HECO + W.POSE + W.PSSP+
  logarea:W.ARTR +logarea:W.HECO + logarea:W.POSE + logarea:W.PSSP+
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("binomial"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

m2c <- inla(survives ~ logarea+ Treatment3 + logarea:Treatment3 + W.ARTR + W.HECO + W.POSE + W.PSSP+
  logarea:W.ARTR +logarea:W.HECO + logarea:W.POSE + logarea:W.PSSP+
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("binomial"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

# compare DIC
tmp <- c("m0","m1a","m1b","m1c","m2a","m2b","m2c") #,"m3a","m3b","m3c")
myDIC <- c(summary(m0)$dic[1],summary(m1a)$dic[1],summary(m1b)$dic[1],summary(m1c)$dic[1],
           summary(m2a)$dic[1],summary(m2b)$dic[1],summary(m2c)$dic[1]) #,
          # summary(m3a)$dic[1],summary(m3b)$dic[1],summary(m3c)$dic[1])
names(myDIC)<-tmp
myDIC
