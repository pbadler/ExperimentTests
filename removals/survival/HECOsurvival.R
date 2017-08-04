# PBA March 2016

# call from removal_analaysis_wrapper.r

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- "HECO"
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste(root,"/ExperimentTests/data/idaho",sep="")
dataDir2 <- paste(root,"/ExperimentTests/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------

#dists <- read.csv(paste(dataDir2,"/speciesdata/IdahoModDistanceWeights_noExptl.csv",sep=""));
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

# import old data--------------------------------------------------------

source("fetchSurvData.r")

D1 <- fetchSdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"

# import modern data--------------------------------------------------------
D2 <- fetchSdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)

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

# Set up ID variables for INLA random effects
allD$GroupID <- as.numeric(allD$Group)
allD$yearID <- 100+as.numeric(allD$year) # for random year offset on intercept

# Treatment effects
m1 <- inla(survives ~ logarea+ Treatment + W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
#  logarea:W.ARTR +logarea:W.HECO + logarea:W.POSE + logarea:W.PSSP+ logarea:W.allcov + logarea:W.allpts +
  f(yearID, model="iid", prior="normal",param=c(0,0.001))+
  f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
  family=c("binomial"), verbose=FALSE,
  control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
  control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

# #treatment effect not significant, fit simpler model
# m.best <- inla(survives ~ logarea+ W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
#   logarea:W.ARTR +logarea:W.HECO + logarea:W.POSE + logarea:W.PSSP+logarea:W.allcov + logarea:W.allpts +
#   f(yearID, model="iid", prior="normal",param=c(0,0.001))+
#   f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
#   f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
#   family=c("binomial"), verbose=FALSE,
#   control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
#   control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

