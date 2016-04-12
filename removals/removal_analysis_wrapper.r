
# call all scripts for removal experiment analysis

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

# 1. get treatment trends #################################
source("treatment_trends_removals.r")

# 2. fit vital rate regressions ###########################

# table to store Treatment effects
trtTests <- data.frame("species"="c","stage"="c","effect"=1,"CI.02.5"=1,"CI.97.5"=1,stringsAsFactors = F)

# fit growth models
library(INLA)
setwd("growth")
source("write_params.r") # get function to format and output parameters
for(iSpp in c("ARTR","HECO","POSE","PSSP")){
  source(paste0(iSpp,"growth.r"))
  # write parameters
  formatGrowthPars(m0,paste0(iSpp,"_growth_noTrt.csv"))
  formatGrowthPars(m1,paste0(iSpp,"_growth_Trt.csv"))    
  # save treatment test
  irow <- dim(trtTests)[1]
  trtTests[irow+1,] <- NA
  trtTests[irow+1,1:2] <- c(iSpp,"growth")
  tmp <- grep("Treatment",row.names(m1$summary.fixed))
  trtTests[irow+1,3:5] <- m1$summary.fixed[tmp,c("mean","0.025quant","0.975quant")]
}
setwd("..")

# fit survival models
library(INLA)
setwd("survival")
source("write_params.r") # get function to format and output parameters
for(iSpp in c("ARTR","HECO","POSE","PSSP")){
  source(paste0(iSpp,"survival.r"))
  # write parameters
  formatSurvPars(m0,paste0(iSpp,"_surv_noTrt.csv"))
  formatSurvPars(m1,paste0(iSpp,"_surv_Trt.csv"))   
  # save treatment test
  irow <- dim(trtTests)[1]
  trtTests[irow+1,] <- NA
  trtTests[irow+1,1:2] <- c(iSpp,"survival")
  tmp <- grep("Treatment",row.names(m1$summary.fixed))
  trtTests[irow+1,3:5] <- m1$summary.fixed[tmp,c("mean","0.025quant","0.975quant")]
}
setwd("..")

# fit recruitment model
library(boot)
library(R2WinBUGS)
setwd("recruitment")
source("call_recruit_m0.r")
source("call_recruit_m1.r")

#save treatment test data for ARTR
irow <- dim(trtTests)[1]
trtTests[irow+1,] <- NA
trtTests[irow+1,1:2] <- c("ARTR","recruitment")
tmp <- grep("intcpt.trt",row.names(out$summary))
trtTests[irow+1,3:5] <- out$summary[tmp[5],c("mean","2.5%","97.5%")]

#save treatment test data for the three grasses
irow <- dim(trtTests)[1]
trtTests[(irow+1):(irow+3),] <- NA
trtTests[(irow+1):(irow+3),1:2] <- cbind(c("HECO","POSE","PSSP"),rep("recruitment",3))
tmp <- grep("intcpt.trt",row.names(out$summary))
trtTests[(irow+1):(irow+3),3:5] <- out$summary[tmp[2:4],c("mean","2.5%","97.5%")]

setwd("..")

# 3. explore neighborhood composition ###################################
setwd("Wdistrib")
source("exploreSurvivalWs.R")
setwd("..")

# 4. get IBM predictions for quadrat cover ###############################

sppList <-  c("ARTR","HECO","POSE","PSSP")

source("validate/get_W_functions.r")  # get neighbor distance decay functions

# do control plots
quadList <- paste0("Q",c(1:6,19:26))
removeSpp <- NULL
trtEffects <- FALSE  # TRUE means use a model that includes removal treatment effects
for(iQuad in quadList){
  qName=iQuad
  doGroup=1
  source("validate/ibm_validate_removal.r")
  source("validate/ibm_validate_removal_1step.r")
}

# do no grass plots
quadList <- c("Q48","Q49","Q51","Q55","Q57","Q58","Q60","Q62")  # no grass
removeSpp <- c("HECO","POSE","PSSP")
for(iQuad in quadList){
  qName=iQuad
  doGroup=1
  trtEffects <- FALSE  
  source("validate/ibm_validate_removal.r")
  source("validate/ibm_validate_removal_1step.r")
  trtEffects <- TRUE  
  source("validate/ibm_validate_removal.r")
  source("validate/ibm_validate_removal_1step.r")
}

# do no shrub plots
quadList <- c("Q47","Q50","Q52","Q53","Q54","Q56","Q59","Q61")  # no shrub
removeSpp <- c("ARTR")
for(iQuad in quadList){
  qName=iQuad
  doGroup=1
  trtEffects <- FALSE  
  source("validate/ibm_validate_removal.r")
  source("validate/ibm_validate_removal_1step.r")
  trtEffects <- TRUE  
  source("validate/ibm_validate_removal.r")
  source("validate/ibm_validate_removal_1step.r")
}

# make figure for simulation results
source("validate/summarize_validate_sims1step.r")

