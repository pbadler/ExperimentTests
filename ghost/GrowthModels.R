# Some exploratory analyses 
# much of this is stolen from ExperimentTests/

rm(list=ls(all=TRUE)); graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/ghost/",sep="")); # modify as needed 

# required packages
library(lme4)
library(boot)
library(mgcv)

source("fetchGrowthData.R")
source("TrimQuadrats.R"); 

### growth model for PSSP 
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
doSpp <- "PSSP"

dataDir1 <- paste(root,"/ExperimentTests/data/idaho",sep="")
dataDir2 <- paste(root,"/ExperimentTests/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights for neighborhood competition ------------------
dists <- read.csv(paste(dataDir2,"/speciesData/IdahoModDistanceWeights.csv",sep=""));
for(j in 1:4) dists[,j] = dists[,j]/sum(dists[,j]); 
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

# import historical data--------------------------------------------------------
D <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D$Treatment <- "Control"
D <- trimQuadrats(D)$data;  # remove records where age==1 is not reliable 
D<-subset(D,age>1);         # remove seedlings (age=1) because 'seedlings are different' 
D<-subset(D,D$logarea.t0> log(0.251)); # remove 'too small to measure' 

m1.gam <- gam(logarea.t1~logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts + s(year,bs="re"), data=D);
summary(m1.gam); # leave out POSE and random slope 
m1.gam <- gam(logarea.t1~logarea.t0 + W.ARTR + W.HECO + W.PSSP + W.allcov + W.allpts + s(year,bs="re"), data=D);

m2.gam <- gam(logarea.t1~logarea.t0 + W.ARTR + W.HECO + W.PSSP + W.allcov + W.allpts + 
              W.PSSP*W.ARTR + W.PSSP*W.HECO + I(W.PSSP^2) + s(year,bs="re"),data=D) 
summary(m2.gam); # eliminate some terms 

m3.gam <- gam(logarea.t1~logarea.t0 + W.ARTR + W.HECO + W.PSSP + W.allcov + W.allpts + 
              W.PSSP*W.ARTR + s(year,bs="re"),data=D) 

summary(m3.gam); AIC(m1.gam,m3.gam); 

### growth model for POSE 
doSpp <- "POSE"

# import historical data--------------------------------------------------------
D <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D$Treatment <- "Control"
D <- trimQuadrats(D)$data;  # remove records where age==1 is not reliable 
D<-subset(D,age>1);         # remove seedlings (age=1) because 'seedlings are different' 
D<-subset(D,D$logarea.t0> log(0.251)); # remove 'too small to measure' 

m1.gam <- gam(logarea.t1~logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts + s(logarea.t0,year,bs="re") + s(year,bs="re"), data=D);
m2.gam <- gam(logarea.t1~logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts + 
              W.POSE*W.ARTR + W.POSE*W.HECO + W.PSSP*W.POSE + I(W.POSE^2) + s(logarea.t0,year,bs="re") + s(year,bs="re"),data=D) 

AIC(m1.gam,m2.gam); 
summary(m2.gam); 
