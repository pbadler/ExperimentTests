# PBA July 2018

# call from ??.r

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- iSpp
dataDir1 <- paste(root,"/ExperimentTests/data/idaho",sep="")
dataDir2 <- paste(root,"/ExperimentTests/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# import old data--------------------------------------------------------

source("analysis/growth/fetchGrowthData.r")

D1 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"

# import modern data--------------------------------------------------------

D2 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)

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
allD$year <- as.factor(allD$year)

# remove ARTR outliers (large plants that obviously do not turn into tiny plants)
if(doSpp=="ARTR"){
  tmp=which(allD$quad=="Q23" & allD$year==1945 & allD$trackID==67)
  tmp=c(tmp,which(allD$quad=="Q12" & allD$year==1955 & allD$trackID==25))
  tmp=c(tmp,which(allD$quad=="Q26" & allD$year==1945 & allD$trackID==73))
  allD=allD[-tmp,]
}


#########################################
#  2. Fit models
#########################################

# use lmer
#library(lme4)

# treatment effect
m1 <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+  
             (1|Group)+(logarea.t0|year),data=allD)


# # use INLA
# # Set up ID variables for INLA random effects
# allD$GroupID <- as.numeric(allD$Group)
# allD$yearID <- 100+as.numeric(allD$year) # for random year offset on intercept
# 
# # Basic treatment effect (intercept) model
# m1 <- inla(logarea.t1 ~ logarea.t0 + Treatment + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts +
#   f(yearID, model="iid", prior="normal",param=c(0,0.001))+
#   f(GroupID, model="iid", prior="normal",param=c(0,0.001))+
#   f(year, logarea.t0, model="iid", prior="normal",param=c(0,0.001)), data=allD,
#   family=c("gaussian"), verbose=FALSE,
#   control.predictor = list(link = 1),control.compute=list(dic=T,mlik=T),
#   control.inla = list(h = 1e-10),Ntrials=rep(1,nrow(allD)))

