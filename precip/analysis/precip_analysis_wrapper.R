
# call all scripts for precip experiment analysis

rm(list=ls(all=TRUE))
graphics.off();

###
### 1. get treatment trends #################################
###

library(texreg) # to save output
library(lme4)

statsOutput <- paste0(getwd(),"/stats_tables.tex")
#source("treatment_trends_removals.r")

# make climate figure
#source("climate_fig.r")

# clean up
tmp=ls() ; tmp=tmp[tmp!="root" & tmp!="statsOutput"]
rm(list=tmp)

###
### 2. fit vital rate regressions ###########################
###

library(texreg) # to save output
library(xtable)
library(lme4)
library(INLA)

# table to store Treatment effects
trtTests <- data.frame("species"="c","stage"="c","effect"=1,"CI.02.5"=1,"CI.97.5"=1,stringsAsFactors = F)

# read in distance weights
dists <- read.csv("~/driversdata/data/idaho_modern/speciesData/IdahoModDistanceWeights_noExptl.csv")

###
### fit survival models (takes ~ 10 minutes)
###

###
### fit growth models (takes < 5 mins)
###


###
### fit recruitment models (this can take a few hours)
###


###
### write results to file
###


# make treatment effect figure

# clean up
tmp=ls() ; tmp=tmp[tmp!="root" & tmp!="statsOutput"]
rm(list=tmp)

###
### 3. explore neighborhood composition ###################################
###


# clean up
tmp=ls() ; tmp=tmp[tmp!="root" & tmp!="statsOutput"]
rm(list=tmp)

###
### 4. get IBM predictions for quadrat cover ###############################
###

sppList <-  c("ARTR","HECO","POSE","PSSP")

# read in distance weights
dists <- read.csv("~/driversdata/data/idaho_modern/speciesData/IdahoModDistanceWeights_noExptl.csv")

#max.CI <- F  # TRUE means use maximum removal effect

# do contemporary control plots
quadList <- paste0("Q",c(1:6,19:26))  
groupList <- c(rep(1,6),rep(6,4),rep(3,4))
removeSpp <- NULL
trtEffects <- FALSE  # TRUE means use a model that includes removal treatment effects
for(iQuad in 1:length(quadList)){
  qName=quadList[iQuad]
  doGroup=groupList[iQuad]
  source("analysis/ibm/ibm_validate_removal.r")   # project forward from 2011
}

# do no grass plots
quadList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
removeSpp <- c("HECO","POSE","PSSP")
for(iQuad in quadList){
  qName=iQuad
  doGroup=1
  trtEffects <- FALSE  
  source("ibm/ibm_validate_removal.r")
  source("ibm/ibm_validate_removal_1step.r")
  trtEffects <- TRUE  
  source("ibm/ibm_validate_removal.r")
  source("ibm/ibm_validate_removal_1step.r")
}

# do no shrub plots
quadList <- c("Q47","Q50","Q52","Q53","Q54","Q56","Q59","Q61")  # no shrub
removeSpp <- c("ARTR")
for(iQuad in quadList){
  qName=iQuad
  doGroup=1
  trtEffects <- FALSE  
  source("ibm/ibm_validate_removal.r")
  source("ibm/ibm_validate_removal_1step.r")
  trtEffects <- TRUE  
  source("ibm/ibm_validate_removal.r")
  source("ibm/ibm_validate_removal_1step.r")
}

# make figure for simulation results
source("ibm/summarize_validate_sims1step.r")












# clean up
tmp=ls() ; tmp=tmp[tmp!="root" & tmp!="statsOutput"]
rm(list=tmp)

###
### 5. get equilibrium cover from an IPM ###############################
###

max.CI <- F  # TRUE means use maximum removal effect

sppList <-  c("ARTR","HECO","POSE","PSSP")

source("analysis/ipm/get_W_functions.r")  # get neighbor distance decay functions

#no treatment effects, all species
init.species <- c(1:4)
tlimit <- 2500
burn.in <- 500
trtEffects=F
source("analysis/ipm/IPM-setup.r")
source("analysis/ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"analysis/ipm/baselineCover.csv",row.names=F)
meanCover1 <- meanCover

#no climate effects
init.species <- c(2:4)
source("analysis/ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"analysis/ipm/baselineCover-noARTR.csv",row.names=F)
meanCover2 <- meanCover

#climate effects 
init.species <- c(2:4)
trtEffects=T
max.CI=F
source("analysis/ipm/IPM-setup.r")
source("analysis/ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"analysis/ipm/removalCover-noARTR.csv",row.names=F)
meanCover3 <- meanCover

#treatment effects 
init.species <- c(2:4)
trtEffects=T
max.CI=T
source("analysis/ipm/IPM-setup.r")
source("analysis/ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"analysis/ipm/removalCover-noARTR-maxCI.csv",row.names=F)
meanCover4 <- meanCover

# # removal treatment effects, ARTR removal, no PSSP, what happens to HECO and POSE?
# init.species <- c(2:3)
# trtEffects=T
# max.CI=F
# source("ipm/IPM-setup.r")
# source("ipm/IPM-getEquilibrium.r")
# write.csv(covSave[(burn.in+1):tlimit,],"ipm/removalCover-noARTRnoPSSP.csv",row.names=F)

simResults <- rbind(meanCover1,meanCover2,meanCover3) 
colnames(simResults) <- sppList
write.csv(simResults,"analysis/ipm/simResults-meanCover.csv",row.names=F)

simResults <- rbind(meanCover1,meanCover2,meanCover4) 
colnames(simResults) <- sppList
write.csv(simResults,"analysis/ipm/simResults-meanCover-maxCI.csv",row.names=F)

#no treatment effects, all species, boost ARTR cover
init.species <- c(1:4)
tlimit <- 2500
burn.in <- 500
trtEffects=F
max.CI=F
source("analysis/ipm/IPM-setup.r")
Rpars$intcpt.yr[,1] <- Rpars$intcpt.yr[,1]+1
source("analysis/ipm/IPM-getEquilibrium.r")
print(rbind(meanCover1,meanCover)) # compare baseline run with this one

simFile <- "analysis/ipm/simResults-meanCover.csv"

source("analysis/ipm/IPM-figures.r")

