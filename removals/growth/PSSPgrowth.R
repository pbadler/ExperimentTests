# PBA March 2016

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/growth",sep="")); # modify as needed 

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- "PSSP"
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste(root,"/driversdata/data/idaho",sep="")
dataDir2 <- paste(root,"/driversdata/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------

dists <- read.csv(paste(dataDir1,"/speciesdata/IdahoDistanceWeights.csv",sep=""));
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

# simplest model
m0 <- lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 

# put indicators on intercept only
m1a <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 
m1b <- lmer(logarea.t1~logarea.t0+Treatment2+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 
m1c <- lmer(logarea.t1~logarea.t0+Treatment3+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 

# put indicators on intercept and size only
m2a <- lmer(logarea.t1~logarea.t0+Treatment+logarea.t0:Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 
m2b <- lmer(logarea.t1~logarea.t0+Treatment2+logarea.t0:Treatment2+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 
m2c <- lmer(logarea.t1~logarea.t0+Treatment3+logarea.t0:Treatment3+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 

# # put indicators on intercept and W's only
# m3a <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+
#               Treatment:W.ARTR + Treatment:W.HECO + Treatment:W.POSE + Treatment:W.PSSP+
#              (1|Group)+(logarea.t0|year),data=allD) 
# m3b <- lmer(logarea.t1~logarea.t0+Treatment2+W.ARTR + W.HECO + W.POSE + W.PSSP+
#               Treatment2:W.ARTR + Treatment2:W.HECO + Treatment2:W.POSE + Treatment2:W.PSSP+
#              (1|Group)+(logarea.t0|year),data=allD) 
# m3c <- lmer(logarea.t1~logarea.t0+Treatment3+W.ARTR + W.HECO + W.POSE + W.PSSP+
#               Treatment3:W.ARTR + Treatment3:W.HECO + Treatment3:W.POSE + Treatment3:W.PSSP+
#              (1|Group)+(logarea.t0|year),data=allD) 

# compare AIC
tmp <- c("m0","m1a","m1b","m1c","m2a","2b","m2c") #,"m3a","3b","m3c")
myAIC <- c(AIC(m0),AIC(m1a),AIC(m1b),AIC(m1c),AIC(m2a),AIC(m2b),AIC(m2c)) #,AIC(m3a),AIC(m3b),AIC(m3c))
names(myAIC)<-tmp
myAIC

# try another subset
tmpD <- subset(allD,Treatment3=="Control")
m0.old <- lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
             (logarea.t0|year),data=tmpD) 
# all controls
tmpD <- subset(allD,Treatment=="Control")
m0.controls <- lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
                      (logarea.t0|year),data=tmpD,control=lmerControl(optimizer="bobyqa")) 

# does effect diminish with time?
allD$trtYears <- as.factor(ifelse(allD$Treatment=="No_shrub",
                       as.numeric(as.character(allD$year))-2010,0))
test <-lmer(logarea.t1~trtYears+logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+W.allcov + W.allpts+
             (logarea.t0|year),data=allD) 
