
rm(list = ls() ) 

setwd("C:/Repos/ExperimentTests/precip/")
root <- "c:/repos/" # this is needed to get to the driversdata directory

#library(texreg) # to save output
library(xtable)
library(lme4)
library(dplyr)
library(glmnet) # package for statistical regularization

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

### elastic net ------------------------------------------------------------

for(i in 1:length(sppList)){
  
  doSpp <- sppList[i]
  
  ## model Intercept year effects
  # prepare data
  trainD <- yrBetas[[i]] %>% filter(year < 2011, year > 1926) # drop first year because of NAs in covariates
  y <- trainD$Intercept
  X <- trainD[,7:NCOL(trainD)]
  X <- scale(X, center = TRUE, scale = TRUE)
  pen_facts <- rep(1, ncol(X)) # penalize all covariates
  lambdas <- 10^seq(2, -2, by = -.005) # sequence of penalties to test
  
  
  enet_out <- cv.glmnet(x = X, 
                       y = y, 
                       lambda = lambdas,
                       penalty.factor = pen_facts,
                       family = "gaussian", 
                       alpha = 0.5, # 0 for ridge, 1 for lasso 
                       standardize = FALSE, 
                       type.measure = "mse",
                       nfolds = 6)
  
  # Collect results into data frames ---------------------------------------------
  cv_scores <- enet_out$cvm
  all_coefs <- as.data.frame(as.matrix(t(enet_out$glmnet.fit$beta)[,1:ncol(X)]))
  colnames(all_coefs) <- colnames(X)
  all_coefs <- all_coefs %>%
  mutate(lambda = log(lambdas)) %>%
  gather(covariate, value, -lambda)

  mse_df <- data.frame(lambda = log(lambdas),
                     score = cv_scores)
  best_lambda <- min(mse_df$lambda[which(mse_df$score == min(mse_df$score))])

  
}

