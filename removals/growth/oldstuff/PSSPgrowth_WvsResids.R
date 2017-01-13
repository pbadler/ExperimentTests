# PBA March 2016

# call from removal_analysis_wrapper.r

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- "PSSP"
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste(root,"/driversdata/data/idaho",sep="")
dataDir2 <- paste(root,"/driversdata/data/idaho_modern",sep="")
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

# look at biggest residuals
m1 <- lmer(logarea.t1~logarea.t0+ Treatment+ W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD)
keep <- which(allD$Treatment=="No_shrub")
check<-allD[keep,]
check$resid <- resid(m1)[keep]
check<-check[order(check$resid,decreasing=TRUE),]

# fit again without PSSP --> Poa mistake
bad <- which(allD$quad=="Q54" & allD$year==2012 & allD$trackID==9)
m1.noOutlier <- lmer(logarea.t1~logarea.t0+ Treatment+ W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD[-bad,])

# explore treatment effects
m0 <- lmer(logarea.t1~logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 

m1 <- lmer(logarea.t1~logarea.t0+ Treatment + W.HECO + W.POSE + W.PSSP+ W.ARTR + W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 

m2 <- lmer(logarea.t1~logarea.t0+ Treatment + W.HECO + W.POSE + W.PSSP+ W.ARTR + W.allcov + W.allpts+
           W.PSSP:Treatment+  
             (1|Group)+(logarea.t0|year),data=allD) 

m3 <- lmer(logarea.t1~logarea.t0+ Treatment + W.HECO + W.POSE + W.PSSP+ W.ARTR + W.allcov + W.allpts+
           W.POSE:Treatment+ W.HECO:Treatment+W.PSSP:Treatment+ 
             (1|Group)+(logarea.t0|year),data=allD) 
print(c(AIC(m0),AIC(m1),AIC(m2),AIC(m3)))

# look at residuals vs marginal PSSP effects
m0 <- lmer(logarea.t1~logarea.t0+ W.ARTR+ W.HECO + W.POSE + W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD) 
png("PSSP_marginalWPSSP.png",height=3.5,width=5,units="in",res=400)
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,3,1,1))
  plot(allD$W.PSSP[allD$Treatment=="Control"],residuals(m0)[allD$Treatment=="Control"],
        xlab="W.PSSP",ylab="Residuals",col="darkgrey")
  abline(lm(residuals(m0)[allD$Treatment=="Control"]~allD$W.PSSP[allD$Treatment=="Control"]),col="darkgrey",lwd=2)
  abline(h=0,lty="dashed",col="black")
  points(allD$W.PSSP[allD$Treatment=="No_shrub"],residuals(m0)[allD$Treatment=="No_shrub"],col="red")
  abline(lm(residuals(m0)[allD$Treatment=="No_shrub"]~allD$W.PSSP[allD$Treatment=="No_shrub"]),col="red",lwd=2)
dev.off()

# look at residuals vs marginal W.ARTR effects
m0 <- lmer(logarea.t1~logarea.t0+ W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts+
             (1|Group)+(logarea.t0|year),data=allD)
png("PSSP_marginalWARTR.png",height=3.5,width=5,units="in",res=400)
  par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,3,1,1))
  plot(allD$W.ARTR[allD$Treatment=="Control"],residuals(m0)[allD$Treatment=="Control"],
       xlab="W.ARTR",ylab="Residuals",col="darkgrey")
  abline(lm(residuals(m0)[allD$Treatment=="Control"]~allD$W.ARTR[allD$Treatment=="Control"]),col="black",lwd=2)
  abline(h=0,lty="dashed",col="black")
  points(allD$W.ARTR[allD$Treatment=="No_shrub"],residuals(m0)[allD$Treatment=="No_shrub"],col="red")
  points(0,mean(residuals(m0)[allD$Treatment=="No_shrub"]),pch=16,col="blue")
dev.off()

#boxplots of residuals by treatment at W.ARTR=0
allD$resids <- residuals(m0)
boxplot(resids~Treatment,data=allD,subset=W.ARTR==0)
summary(lm(resids~Treatment,data=allD,subset=W.ARTR==0))
