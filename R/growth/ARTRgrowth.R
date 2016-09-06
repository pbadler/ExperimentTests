# PBA March 2016
# modified by ARK July 2016
# call from precip_analysis_wrapper.r
# simple test for differences in growth rate between precip treatments in the contemporary period 

rm(list = ls())

#########################################
#  1. Import data and calculate W's
#########################################

doSpp <- "ARTR"
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste("~/driversdata/data/idaho",sep="")
dataDir2 <- paste("~/driversdata/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------
#dists <- read.csv('~/driversdata/data/idaho/speciesData/IdahoDistanceWeights.csv')
dists <- read.csv('~/driversdata/data/idaho_modern/speciesData/IdahoModDistanceWeights_noExptl.csv')
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

# import old data--------------------------------------------------------

source("R/growth/fetchGrowthData.r")

D1 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir1,distWts=dists)
D1$Treatment <- "Control"
D1$Period <- "Historical"

# import modern data--------------------------------------------------------
D2 <- fetchGdat(doSpp=doSpp,speciesList=sppList,datadir=dataDir2,distWts=dists)
D2$Period <- "Modern"

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
  keep <- which(is.element(allD$Treatment,c("Control","No_grass", "Irrigation", "Drought")))
}else{
  keep <- which(is.element(allD$Treatment,c("Control","No_shrub", "Irrigation", "Drought")))
}
allD <- allD[keep,]

# remove outliers (large plants that obviously do not turn into tiny plants)
tmp=which(allD$quad=="Q23" & allD$year==1945 & allD$trackID==67)
tmp=c(tmp,which(allD$quad=="Q12" & allD$year==1955 & allD$trackID==25))
tmp=c(tmp,which(allD$quad=="Q26" & allD$year==1945 & allD$trackID==73))
allD=allD[-tmp,]

#########################################
#  2. Fit models
#########################################

# set up indicator variables
allD$Treatment2 <- allD$Treatment
allD$Treatment2[allD$Treatment=="Control" & allD$year>2000] <- "ControlModern"
allD$Treatment[ allD$year < 2012 & allD$Treatment %in% c('Drought', 'Irrigation') ] <- 'Control'  # set initial treatment to control

allD <- subset(allD, Treatment %in% c('Control', 'Drought', 'Irrigation') )

allD$year <- as.factor(allD$year)

# use lmer
library(lme4)
# w/o treatment effect

m0 <- lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
              (1|Group)+(logarea.t0|year),data=subset(allD, as.numeric(levels(year)[year]) > 2006)) 

# w/ treatment effect
m1 <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
              (1|Group)+(logarea.t0|year),data=subset(allD, as.numeric(levels(year)[year]) > 2006)) 
# 
anova(m1, m0) # no treatment effect 

lmer_results = list(m0, m1)

saveRDS(lmer_results, file = 'output/ARTR_growth_treatment_effects.lmer.RDS')
