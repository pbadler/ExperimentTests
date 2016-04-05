
# call all scripts for removal analysis

rm(list=ls(all=TRUE))
graphics.off();

root=ifelse(.Platform$OS.type=="windows","c:/Repos","~/repos"); # modify as needed
setwd(paste(root,"/ExperimentTests/removals/",sep="")); # modify as needed 

# get treatment trends
source("treatment_trends_removals.r")

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
library("INLA")
setwd("survival")
#source("write_params.r") # get function to format and output parameters
for(iSpp in c("ARTR","HECO","POSE","PSSP")){
  source(paste0(iSpp,"survival.r"))
  #formatGrowthPars(m0,paste0(iSpp,"_growth_noTrt.csv"))
  #formatGrowthPars(m1,paste0(iSpp,"_growth_Trt.csv"))       
}
setwd("..")
