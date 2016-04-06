
# call all scripts for removal experiment analysis

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

# 1. get treatment trends #################################
source("treatment_trends_removals.r")

# 2. fit vital rate regressions ###########################

# fit growth models
library(lme4)
setwd("growth")
source("write_params.r") # get function to format and output parameters
for(iSpp in c("ARTR","HECO","POSE","PSSP")){
  source(paste0(iSpp,"growth.r"))
  formatGrowthPars(m0,paste0(iSpp,"_growth_noTrt.csv"))
  formatGrowthPars(m1,paste0(iSpp,"_growth_Trt.csv"))       
}
setwd("..")

# fit survival models
library(INLA)
setwd("survival")
source("write_params.r") # get function to format and output parameters
for(iSpp in c("ARTR","HECO","POSE","PSSP")){
  source(paste0(iSpp,"survival.r"))
  formatSurvPars(m0,paste0(iSpp,"_surv_noTrt.csv"))
  formatSurvPars(m1,paste0(iSpp,"_surv_Trt.csv"))       
}
setwd("..")

# fit recruitment model
library(boot)
library(R2WinBUGS)
setwd("recruitment")
source("call_recruit_m0.r")
source("call_recruit_m1.r")
setwd("..")

# 3. explore neighborhood composition #################################
setwd("Wdistrib")
source("exploreSurvivalWs.R")
setwd("..")

# get IBM predictions for quadrat cover ###############################



