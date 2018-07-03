
rm(list = ls() ) 

setwd("C:/Repos/ExperimentTests/precip/")
root <- "c:/repos/" # this is needed to get to the driversdata directory

#library(texreg) # to save output
library(xtable)
library(lme4)

# read in distance weights
dists <- read.csv(paste0(root,"/ExperimentTests/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))

# object to save year effects
sppList <- c("ARTR","HECO","POSE","PSSP")
yrBetas <- vector("list",length(sppList))

for(i in 1:length(sppList)){
  iSpp <- sppList[i]
  source("analysis/growth/inla_growth.r")
  yrBetas[[i]] <- ranef(m1)$year
}

