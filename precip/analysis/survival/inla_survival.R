# PBA March 2016

# call from removal_analaysis_wrapper.r

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- iSpp
dataDir1 <- paste(root,"/ExperimentTests/data/idaho",sep="")
dataDir2 <- paste(root,"/ExperimentTests/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# import old data

source("analysis/survival/fetchSurvData.r")

D1 <- fetchSdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"

# import modern data
D2 <- fetchSdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)

# merge in treatment data
tmp <- read.csv(paste(dataDir2,"/quad_info.csv",sep=""))
tmp <- tmp[,c("quad","Treatment")]
D2 <- merge(D2,tmp, all.x=T)

# combine old and modern
allD <- rbind(D1,D2)
rm(D1,D2,tmp)

# clean up dataset ----------------------------------------------
allD$year[allD$year<2000] <- allD$year[allD$year<2000] + 1900

keep <- which(is.element(allD$Treatment,c("Control","Drought","Irrigation")))
allD <- allD[keep,]

#########################################
#  2. Fit models
#########################################
library(INLA)

allD$year <- as.factor(allD$year)

# Set up ID variables for INLA random effects
allD$GroupID <- as.numeric(allD$Group)
allD$yearID <- 100+as.numeric(allD$year) # for random year offset on intercept

m1 <- inla(survives ~ logarea+ Treatment + W.ARTR + W.HECO + W.POSE + W.PSSP +
   f(yearID, model="iid", prior="normal",param=c(0,0.001))+
   f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
   f(year, logarea, model="iid", prior="normal",param=c(0,0.001)), data=allD,
   family=c("binomial"), verbose=FALSE,
   control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
   control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

