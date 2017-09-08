# PBA March 2016

# call from removal_analysis_wrapper.r

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- "POSE"
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste(root,"/ExperimentTests/data/idaho",sep="")
dataDir2 <- paste(root,"/ExperimentTests/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------

#dists <- read.csv(paste(dataDir2,"/speciesdata/IdahoModDistanceWeights_noExptl.csv",sep=""));
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

# import old data--------------------------------------------------------

source("fetchGrowthData.r")

D1 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"

# import modern data--------------------------------------------------------

D2 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)

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
rm(D1,D2,tmp)

# merge data on removals at individual level
tmp <- read.csv(paste(dataDir2,"/speciesData/",doSpp,"/",doSpp,"_within_ARTRremovals.csv",sep=""))
tmp <- tmp[,c("quad","year","trackID","inARTR")] 
allD<-merge(allD,tmp,all.x=T)
allD$inARTR[is.na(allD$inARTR)] <- 0

# clean up dataset ----------------------------------------------
allD$year[allD$year<2000] <- allD$year[allD$year<2000] + 1900

if(doSpp=="ARTR"){
  keep <- which(is.element(allD$Treatment,c("Control","No_grass")))
}else{
  keep <- which(is.element(allD$Treatment,c("Control","No_shrub")))
}
allD <- allD[keep,]

# remove outliers (large plants that obviously do not turn into tiny plants)


#########################################
#  2. Fit models
#########################################

library(lme4)

# set up indicator variables
allD$Treatment2 <- allD$Treatment
allD$Treatment2[allD$year>2000] <- "Modern"
allD$Treatment3 <- allD$Treatment
allD$Treatment3[allD$Treatment=="Control" & allD$year>2000] <- "ControlModern"

allD$year <- as.factor(allD$year)

# use INLA
# Set up ID variables for INLA random effects
allD$GroupID <- as.numeric(allD$Group)
allD$yearID <- 100+as.numeric(allD$year) # for random year offset on intercept

# Treatment effect
m1 <- inla(logarea.t1 ~ logarea.t0 + Treatment + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts +
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
#  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("gaussian"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

# explore additional models

# add individual level removal info to best model
m2 <- inla(logarea.t1 ~ logarea.t0 + Treatment + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts + inARTR +
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
#  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("gaussian"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

# save fixed effects summary to file
cat("",file=statsOutput,sep="\n",append=T)
cat(capture.output(print(xtable(m2$summary.fixed,digits=4,caption=paste("Summary of fixed effects for the Po. secunda model with individual-level A. tripartita removal data"),
        label=paste0(iSpp,"growth-inARTR")),caption.placement="top")),file=statsOutput,sep="\n",append=T)
rm(m2)

# does effect diminish with time?
allD$trtYears <- as.factor(ifelse(allD$Treatment=="No_shrub",
                       as.numeric(as.character(allD$year))-2010,0))
m3 <- inla(logarea.t1 ~ trtYears + logarea.t0  + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts + inARTR +
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
#  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("gaussian"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

# save fixed effects summary to file
cat("",file=statsOutput,sep="\n",append=T)
cat(capture.output(print(xtable(m3$summary.fixed,digits=4,caption=paste("Summary of fixed effects for the Po. secunda model with treatment*year effects"),
        label=paste0(iSpp,"growth-trtYears")),caption.placement="top")),file=statsOutput,sep="\n",append=T)
rm(m3)

# see if parameters change if we only use data from control plots
controlD <- subset(allD,Treatment=="Control")
controlD$GroupID <- as.numeric(controlD$Group)
controlD$yearID <- 100+as.numeric(controlD$year) # for random year offset on intercept
m1.control <- inla(logarea.t1 ~ logarea.t0  + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts +
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001)), data=controlD,
  family=c("gaussian"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(controlD)))
plot(m1.control$summary.fixed$mean,m1$summary.fixed$mean[-3],xlab="Fixed effects (controls)",
     ylab="Fixed effects (all plots)")
abline(0,1)
cor(m1.control$summary.fixed$mean,m1$summary.fixed$mean[-3])

