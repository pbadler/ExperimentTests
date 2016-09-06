# PBA March 2016
# modified by ARK July 2016
# call from removal_analysis_wrapper.r

#########################################
#  1. Import data and calculate W's
#########################################

root <- '~'
doSpp <- "PSSP"
sppList <- c("ARTR","HECO","POSE","PSSP","allcov","allpts")
dataDir1 <- paste(root,"/driversdata/data/idaho",sep="")
dataDir2 <- paste(root,"/driversdata/data/idaho_modern",sep="")
nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 

# set up distance weights------------------------------------------------

dists <- read.csv(paste(dataDir2,"/speciesData/IdahoModDistanceWeights_noExptl.csv",sep=""));
dists$allcov <- rowMeans(dists[,1:4])  # for "other" polygons use average of big 4
dists$allpts <- dists$POSE  # set forb dist wts = smallest grass (POSE)

# import old data--------------------------------------------------------
source("./R/growth/fetchGrowthData.r")

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
  keep <- which(is.element(allD$Treatment,c("Control","No_grass", "Irrigation", "Drought")))
}else{
  keep <- which(is.element(allD$Treatment,c("Control","No_shrub", "Irrigation", "Drought")))
}
allD <- allD[keep,]

# remove outliers (large plants that obviously do not turn into tiny plants)




# set up indicator variables 
allD$Treatment2 <- allD$Treatment
allD$Treatment2[allD$year>2000] <- "Modern"
allD$Treatment3 <- allD$Treatment
allD$Treatment3[allD$Treatment=="Control" & allD$year>2000] <- "ControlModern"
allD$Treatment[ allD$year < 2012 & allD$Treatment %in% c('Drought', 'Irrigation') ] <- 'Control'  # set initial treatment to control

# ----------- use this data for prediction ------------------------------------------------------------------------------

saveRDS(allD, 'data/temp_data/PSSP_growth.RDS') 

# -----------------------------------------------------------------------------------------------------------------------

#########################################
#  2. Fit models
#########################################

library(lme4)
library(INLA)

allD <- subset(allD, Treatment %in% c('Control', 'Drought', 'Irrigation') )

allD$year <- as.factor(allD$year)
# 
# 

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

saveRDS(lmer_results, file = 'output/PSSP_growth_treatment_effects.lmer.RDS')


#
# output<-capture.output(texreg(m2.lmer, ci.force=TRUE,label="table:PSSPgrowth-inARTR",
#       caption="\textit{P. spicata} growth with \textit{Artemisia} canopy effect",
#       caption.above=TRUE))
# cat(output,file=statsOutput,sep="\n",append=T)
# cat("",file=statsOutput,sep="\n",append=T)
# 
# # does effect diminish with time?
# allD$trtYears <- as.factor(ifelse(allD$Treatment=="No_shrub",
#                        as.numeric(as.character(allD$year))-2010,0))
# m1.time <-lmer(logarea.t1~trtYears+logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+ W.allcov + W.allpts +
#              (logarea.t0|year),data=allD) 
# output<-capture.output(texreg(m1.time, ci.force=TRUE,label="table:PSSPgrowth-byYr",
#       caption="\textit{P. spicata} growth with year-by-treatment interaction",
#       caption.above=TRUE))
# 
# cat(output,file=statsOutput,sep="\n",append=T)
# cat("",file=statsOutput,sep="\n",append=T)
# 
