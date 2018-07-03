
rm(list = ls() ) 

setwd("C:/Repos/ExperimentTests/precip/")
root <- "c:/repos/" # this is needed to get to the driversdata directory

#library(texreg) # to save output
library(xtable)
library(lme4)

# read in distance weights
dists <- read.csv(paste0(root,"/ExperimentTests/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))

# import climate covariates
Cdat <- readRDS('data/temp_data/all_clim_covs.RDS')

# object to save year effects
sppList <- c("ARTR","HECO","POSE","PSSP")
yrBetas <- vector("list",length(sppList))

for(i in 1:length(sppList)){
 
   iSpp <- sppList[i]
  # fit demography model
  source("analysis/growth/inla_growth.r")
   
  #extract parameters for controls
  yrBetas[[i]] <- ranef(m1)$year
  yrBetas[[i]]$year <- row.names(yrBetas[[i]])
  yrBetas[[i]]$Treatment <- "Control"
  names(yrBetas[[i]])[1] <- "Intercept"
  
  #add parameters for drought treatment
  tmp <- subset(yrBetas[[i]],year > 2010)
  tmp$Treatment <- "Drought"
  tmp$Intercept <- tmp$Intercept + fixef(m1)[which(names(fixef(m1))=="TreatmentDrought")]
  yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
  
  #add parameters for irrigation treatment
  tmp <- subset(yrBetas[[i]],year > 2010 & Treatment=="Control")
  tmp$Treatment <- "Irrigation"
  tmp$Intercept <- tmp$Intercept + fixef(m1)[which(names(fixef(m1))=="TreatmentIrrigation")]
  yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
  
  #merge in climate covariates
  yrBetas[[i]] <- merge(yrBetas[[i]],Cdat,all.x=T)
  
}



