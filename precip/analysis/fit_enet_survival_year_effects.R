
rm(list = ls() ) 

setwd("C:/Repos/ExperimentTests/precip/")
root <- "c:/repos/" # this is needed to get to the driversdata directory

#library(texreg) # to save output
library(xtable)
library(lme4)
library(dplyr)
library(glmnet) # package for statistical regularization
library(INLA)
library(boot)
source("analysis/figure_scripts/elastic_net_observed_vs_predicted.R")

# ### 1. Fit models ------------------------------------------------------
# 
# # read in distance weights
# dists <- read.csv(paste0(root,"/ExperimentTests/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))
# 
# # import climate covariates
# Cdat <- readRDS('data/temp_data/all_clim_covs.RDS')
# 
# # object to save year effects
# sppList <- c("ARTR","HECO","POSE","PSSP")
# yrBetas <- vector("list",length(sppList))
# 
# for(i in 1:length(sppList)){
#  
#    iSpp <- sppList[i]
#   # fit demography model
#   source("analysis/survival/inla_survival.r")
#    
#   #extract parameters for controls
#   Intercept=m1$summary.random$yearID$mean
#   year=m1$summary.random$year$ID
#   logarea=m1$summary.random$year$mean
#   yrBetas[[i]]=data.frame(year,Intercept,logarea)
#   yrBetas[[i]]$Treatment <- "Control"
#   yrBetas[[i]]$year <-  as.numeric(as.character(yrBetas[[i]]$year))
# 
#   #add parameters for drought treatment
#   tmp <- subset(yrBetas[[i]],year > 2010 & Treatment == "Control")
#   tmp$Treatment <- "Drought"
#   tmp$Intercept <- tmp$Intercept + m1$summary.fixed$mean[which(m1$names.fixed=="TreatmentDrought")]
#   yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
#   
#   #add parameters for irrigation treatment
#   tmp <- subset(yrBetas[[i]], year > 2010 & Treatment == "Control")
#   tmp$Treatment <- "Irrigation"
#   tmp$Intercept <- tmp$Intercept + m1$summary.fixed$mean[which(m1$names.fixed=="TreatmentIrrigation")]
#   yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
#   
#   #merge in climate covariates
#   yrBetas[[i]] <- merge(yrBetas[[i]],Cdat,all.x=T)
#   
# }

### 1. import year effects from stan model ------------------------

# import climate covariates
Cdat <- readRDS('data/temp_data/all_clim_covs.RDS')

# # limit to a few climate variables
# keep <- grep("VWC.sp",names(Cdat))
# keep <- c(keep,grep("T.sp",names(Cdat)))
# Cdat <-Cdat[,c(1:4,keep)]

# object to save year effects
sppList <- c("ARTR","HECO","POSE","PSSP")
yrBetas <- vector("list",length(sppList))
size_info <- vector("list",length(sppList))

for(i in 1:length(sppList)){
 
  iSpp <- sppList[i]
  
  # get vector of years
  dat <- readRDS("data/temp_data/survival_data_lists_for_stan.RDS")[[i]]
  year <- unique(dat$year2)
  size_small <-quantile(dat$X,0.1) ; size_large <- quantile(dat$X,0.9) # get 10% and 90% size quantiles
  size_info[[i]] <- list( size_small, size_large, dat$Xcenter, dat$Xscale )
  names(size_info[[i]]) <- c("small","large","center","scale")
  rm(dat)
  
  # import model summary
  m1 <- read.csv(paste0("output/treatment_model_parameters_",iSpp,"_survival.csv"))
   
  #extract parameters for controls
  keep <- grep("a", substring(as.character(m1$X),1,1))
  Intercept <- m1$mean[keep]
  keep <- grep("b1", substring(as.character(m1$X),1,2))
  keep <- keep[m1$X[keep]!="b1_mu"] 
  logarea <- m1$mean[keep]
  yrBetas[[i]] <- data.frame(year,Intercept,logarea)
  yrBetas[[i]]$Treatment <- "Control"
  
  #add parameters for drought treatment
  tmp <- subset(yrBetas[[i]],year > 2010 & Treatment == "Control")
  tmp$Treatment <- "Drought"
  tmp$Intercept <- tmp$Intercept + m1$mean[m1$X=="bt[1]"]
  yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
  
  #add parameters for irrigation treatment
  tmp <- subset(yrBetas[[i]], year > 2010 & Treatment == "Control")
  tmp$Treatment <- "Irrigation"
  tmp$Intercept <- tmp$Intercept + + m1$mean[m1$X=="bt[2]"]
  yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
  
  # calculate survival of small and large plants 
  yrBetas[[i]]$surv_small <- yrBetas[[i]]$Intercept + yrBetas[[i]]$logarea*size_info[[i]]$small
  yrBetas[[i]]$surv_big <- yrBetas[[i]]$Intercept + yrBetas[[i]]$logarea*size_info[[i]]$large
  yrBetas[[i]]$surv_small_resp <- inv.logit(yrBetas[[i]]$surv_small)
  yrBetas[[i]]$surv_big_resp <- inv.logit(yrBetas[[i]]$surv_big)
  
  #merge in climate covariates
  yrBetas[[i]] <- merge(yrBetas[[i]],Cdat,all.x=T)
  
}


### 2. elastic net on intercept year effects ------------------------------------------------------------

enet_predictions <- vector("list",length(sppList))
best_coefs_small <- matrix(NA,length(sppList),ncol(Cdat)-4) # warning: hardwired "4" on climate columns
colnames(best_coefs_small) <- names(Cdat[,5:ncol(Cdat)])
best_coefs_big <- best_coefs_small

for(i in 1:length(sppList)){
  
  doSpp <- sppList[i]
  
  # prepare training data
  trainD <- yrBetas[[i]] %>% filter(year < 2011, year > 1926) # drop first year because of NAs in covariates
  y_small <- trainD$surv_small
  y_big <- trainD$surv_big
  X <- trainD[,11:NCOL(trainD)]
  X_mean <- colMeans(X); X_sd <- apply(X,MARGIN=2,FUN="sd")
  X <- scale(X, center = X_mean, scale = X_sd)
  
  # prepare out of sample data
  newD <- yrBetas[[i]] %>% filter(year >= 2011)
  y_new_small <- newD$surv_small
  y_new_big <- newD$surv_big
  X_new <- newD[,11:NCOL(newD)]
  X_new <- scale(X_new, center = X_mean, scale = X_sd)
  
  # fit elastic net on intercept
  pen_facts <- rep(1, ncol(X)) # penalize all covariates
  lambdas <- 10^seq(2, -2, by = -0.005) # sequence of penalties to test
  enet_small <- cv.glmnet(x = X, 
                       y = y_small, 
                       lambda = lambdas,
                       penalty.factor = pen_facts,
                       family = "gaussian", 
                       alpha = 0.5, # 0 for ridge, 1 for lasso 
                       standardize = FALSE, 
                       type.measure = "mse",
                       nfolds = length(y_small))
  
  best_coefs_small[i,] <- enet_small$glmnet.fit$beta[,which(enet_small$lambda==enet_small$lambda.min)]
  
  # predictions for training data
  y_hat_small <- predict(enet_small,newx=X,s="lambda.min")
  
  # predictions for out of sample data
  y_hat_new_small <- predict(enet_small,newx=X_new, s="lambda.min")

  # fit elastic net on slope year effects
  enet_big <- cv.glmnet(x = X, 
                       y = y_big, 
                       lambda = lambdas,
                       penalty.factor = pen_facts,
                       family = "gaussian", 
                       alpha = 0.5, # 0 for ridge, 1 for lasso 
                       standardize = FALSE, 
                       type.measure = "mse",
                       nfolds = length(y_big))
  
  best_coefs_big[i,] <- enet_big$glmnet.fit$beta[,which(enet_big$lambda==enet_big$lambda.min)]
  
  # predictions for training data
  y_hat_big <- predict(enet_big,newx=X,s="lambda.min")
  
  # predictions for out of sample data
  y_hat_new_big <- predict(enet_big,newx=X_new, s="lambda.min")
  
  # save output
  enet_predictions[[i]] <- list(trts = newD$Treatment,
                                y_small=y_small,y_hat_small=y_hat_small,y_new_small=y_new_small, y_hat_new_small=y_hat_new_small,
                                y_big=y_big, y_hat_big=y_hat_big, y_new_big=y_new_big, y_hat_new_big=y_hat_new_big)
  saveRDS(best_coefs_small,"analysis/survival/enet_best_coefs_small.RDS")
  saveRDS(best_coefs_big,"analysis/survival/enet_best_coefs_big.RDS")

}

# clean up
rm(list= ls()[!(ls() %in% c('best_coefs_big','best_coefs_small','enet_predictions',
                            'sppList','yrBetas','size_info','obs_pred_fig'))]) 


### 3. Figures -----------------------------------------------------------------

# loop through species and plot observed and predicted intercepts and slopes

pdf("figures/enet_survival_yr_effects.pdf",height=4,width=8)
par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1))

for(i in 1:length(sppList)){
  
  # do small plants
  trts <- enet_predictions[[i]]$trts
  keep <- grep("small",names(enet_predictions[[i]]))
  fig_dat <- enet_predictions[[i]][keep]
  names(fig_dat) <- sub("_small","",names(fig_dat)) # trim names
  fig_dat <- lapply(fig_dat,FUN="inv.logit") 
  obs_pred_fig(fig_dat, trts, vital_rate="survival", effect="small plants", legend_location= "bottomright")
  
  # do big plants
  trts <- enet_predictions[[i]]$trts
  keep <- grep("big",names(enet_predictions[[i]]))
  fig_dat <- enet_predictions[[i]][keep]
  names(fig_dat) <- sub("_big","",names(fig_dat)) # trim names
  fig_dat <- lapply(fig_dat,FUN="inv.logit") 
  obs_pred_fig(fig_dat, trts, vital_rate="survival", effect="big plants", legend_location= "bottomright")
  
}

dev.off()
