
# call all scripts for removal experiment analysis

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

###
### 1. get treatment trends #################################
###

library(texreg) # to save output
library(lme4)

statsOutput <- paste0(getwd(),"/stats_tables.tex")
source("treatment_trends_removals.r")

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
dists <- read.csv(paste0(root,"/driversdata/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))

# fit survival models (takes ~ 10 minutes)
setwd("survival")
source("write_params.r") # get function to format and output parameters
for(iSpp in c("ARTR","HECO","POSE","PSSP")){
  
  source(paste0(iSpp,"survival.r"))
  
  # save fixed effects summary to file
  cat("",file=statsOutput,sep="\n",append=T)
  cat(capture.output(print(xtable(m1$summary.fixed,digits=4,caption=paste("Summary of fixed effects for the",iSpp,"survival model"),
        label=paste0(iSpp,"survival")),caption.placement="top")),file=statsOutput,sep="\n",append=T)
  
  # save treatment test
  irow <- dim(trtTests)[1]
  trtTests[irow+1,] <- NA
  trtTests[irow+1,1:2] <- c(iSpp,"survival")
  tmp <- grep("Treatment",row.names(m1$summary.fixed))
  trtTests[irow+1,3:5] <- m1$summary.fixed[tmp,c("mean","0.025quant","0.975quant")]
  
  # write parameters
  formatSurvPars(m1,paste0(iSpp,"_surv.csv")) 
  
}
setwd("..")

# fit growth models (takes < 5 mins)
setwd("growth")
source("write_params.r") # get function to format and output parameters

for(iSpp in c("ARTR","HECO","POSE","PSSP")){
  
  source(paste0(iSpp,"growth.r"))
  
  # save fixed effects summary to file
  cat("",file=statsOutput,sep="\n",append=T)
  cat(capture.output(print(xtable(m1$summary.fixed,digits=4,caption=paste("Summary of fixed effects for the",iSpp,"growth model"),
        label=paste0(iSpp,"growth")),caption.placement="top")),file=statsOutput,sep="\n",append=T)
  
  # save treatment test
  irow <- dim(trtTests)[1]
  trtTests[irow+1,] <- NA
  trtTests[irow+1,1:2] <- c(iSpp,"growth")
  tmp <- grep("Treatment",row.names(m1$summary.fixed))
  trtTests[irow+1,3:5] <- m1$summary.fixed[tmp,c("mean","0.025quant","0.975quant")]
  
  # write parameters for best model
  formatGrowthPars(m1,paste0(iSpp,"_growth.csv")) 
  
}
setwd("..")

# fit recruitment models (this can take a few hours)
library(boot)
library(R2WinBUGS)
setwd("recruitment")
#source("call_recruit_m1.r")

# add treatment test data for ARTR 
pars.summary <- read.csv("recruit_params_m1.csv")
irow <- dim(trtTests)[1]
trtTests[irow+1,] <- NA
trtTests[irow+1,1:2] <- c("ARTR","recruitment")
tmp <- grep("intcpt.trt",row.names(pars.summary))
trtTests[irow+1,3:5] <- pars.summary[tmp[5],c("mean","X2.5.","X97.5.")]

# add treatment test data for the three grasses
irow <- dim(trtTests)[1]
trtTests[(irow+1):(irow+3),] <- NA
trtTests[(irow+1):(irow+3),1:2] <- cbind(c("HECO","POSE","PSSP"),rep("recruitment",3))
trtTests[(irow+1):(irow+3),3:5] <- pars.summary[tmp[2:4],c("mean","X2.5.","X97.5.")]

# create summary stats table
fixeffs <- c("intcpt.mu","intcpt.trt","dd","theta","u")
keep <- NULL
for(i in fixeffs){
  keep <- c(keep,grep(i,row.names(pars.summary),fixed=TRUE))
}
keep <- keep[-c(5,10:12,33:44)]  # some editing
output <- pars.summary[keep,c(1:3,7:9)]
cat("",file=statsOutput,sep="\n",append=T)
cat(capture.output(print(xtable(output,digits=4,caption=paste("Summary of fixed effects for the recruitment model"),
        label="table:recruitment"),caption.placement="top")),file=statsOutput,sep="\n",append=T)
  
setwd("..")

# write trtTests to file
trtTests <- trtTests[-1,] # throw away first line junk
write.csv(trtTests,"treatment_test_results.csv",row.names=F)

# make treatment effect figure
source("treatment_test_figure.r")

###
### 3. explore neighborhood composition ###################################
###
library("TeachingDemos") # for inset plots

source("exploreSurvivalWs.R")


###
### 4. get IBM predictions for quadrat cover ###############################
###

sppList <-  c("ARTR","HECO","POSE","PSSP")

# read in distance weights
dists <- read.csv(paste0(root,"/driversdata/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))




source("ibm/get_W_functions.r")  # get neighbor distance decay functions

# do contemporary control plots
quadList <- paste0("Q",c(1:6,19:26))  
groupList <- c(rep(1,6),rep(6,4),rep(3,4))
removeSpp <- NULL
trtEffects <- FALSE  # TRUE means use a model that includes removal treatment effects
for(iQuad in 1:length(quadList)){
  qName=quadList[iQuad]
  doGroup=groupList[iQuad]
  source("ibm/ibm_validate_removal.r")   # project forward from 2011
  source("ibm/ibm_validate_removal_1step.r")  # just predict one time step ahead for each year
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


###
### 5. get equilibrium cover from an IPM ###############################
###

sppList <-  c("ARTR","HECO","POSE","PSSP")

source("ibm/get_W_functions.r")  # get neighbor distance decay functions

#no treatment effects, all species
init.species <- c(1:4)
tlimit <- 2500
burn.in <- 500
trtEffects=F
source("ipm/IPM-setup.r")
source("ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"ipm/baselineCover.csv",row.names=F)
meanCover1 <- meanCover

#no treatment effects, ARTR removal
init.species <- c(2:4)
source("ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"ipm/baselineCover-noARTR.csv",row.names=F)
meanCover2 <- meanCover

# removal treatment effects, ARTR removal
init.species <- c(2:4)
trtEffects=T
source("ipm/IPM-setup.r")
source("ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"ipm/removalCover-noARTR.csv",row.names=F)
meanCover3 <- meanCover

# removal treatment effects, ARTR removal, no PSSP
init.species <- c(2:3)
trtEffects=T
source("ipm/IPM-setup.r")
source("ipm/IPM-getEquilibrium.r")
write.csv(covSave[(burn.in+1):tlimit,],"ipm/removalCover-noARTRnoPSSP.csv",row.names=F)
meanCover4 <- meanCover

simResults <- rbind(meanCover1,meanCover2,meanCover3,meanCover4)
colnames(simResults) <- sppList
write.csv(simResults,"ipm/simResults-meanCover.csv",row.names=F)

source("ipm/IPM-figures.r")

