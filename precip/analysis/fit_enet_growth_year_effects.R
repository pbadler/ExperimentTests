
rm(list = ls() ) 

setwd("C:/Repos/ExperimentTests/precip/")
root <- "c:/repos/" # this is needed to get to the driversdata directory

#library(texreg) # to save output
library(xtable)
library(lme4)
library(dplyr)
library(glmnet) # package for statistical regularization
source("analysis/figure_scripts/elastic_net_observed_vs_predicted.R")

# ### 1. import year effects from stan model ------------------------
# 
# # import climate covariates
# Cdat <- readRDS('data/temp_data/all_clim_covs.RDS')
# 
# # object to save year effects
# sppList <- c("ARTR","HECO","POSE","PSSP")
# yrBetas <- vector("list",length(sppList))
# size_info <- vector("list",length(sppList))
# 
# for(i in 1:length(sppList)){
# 
#   iSpp <- sppList[i]
# 
#   # get vector of years
#   dat <- readRDS("data/temp_data/growth_data_lists_for_stan.RDS")[[i]]
#   year <- unique(dat$year2)
#   size_small <-quantile(dat$X,0.1) ; size_large <- quantile(dat$X,0.9) # get 10% and 90% size quantiles
#   size_info[[i]] <- list( size_small, size_large, dat$Xcenter, dat$Xscale )
#   names(size_info[[i]]) <- c("small","large","center","scale")
#   rm(dat)
# 
#   # import model summary
#   m1 <- read.csv(paste0("output/treatment_model_parameters_",iSpp,"_growth.csv"))
# 
#   #extract parameters for controls
#   keep <- grep("a", substring(as.character(m1$X),1,1))
#   Intercept <- m1$mean[keep]
#   keep <- grep("b1", substring(as.character(m1$X),1,2))
#   keep <- keep[m1$X[keep]!="b1_mu"] # exclude mean slope
#   logarea.t0 <- m1$mean[keep]
#   yrBetas[[i]] <- data.frame(year,Intercept,logarea.t0)
#   yrBetas[[i]]$Treatment <- "Control"
# 
#   #add parameters for drought treatment
#   tmp <- subset(yrBetas[[i]],year > 2010 & Treatment == "Control")
#   tmp$Treatment <- "Drought"
#   tmp$Intercept <- tmp$Intercept + m1$mean[m1$X=="bt[1]"]
#   yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
# 
#   #add parameters for irrigation treatment
#   tmp <- subset(yrBetas[[i]], year > 2010 & Treatment == "Control")
#   tmp$Treatment <- "Irrigation"
#   tmp$Intercept <- tmp$Intercept + m1$mean[m1$X=="bt[2]"]
#   yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
# 
#   # calculated relative growth rates of small and large plants on arithmetic scale
#   yrBetas[[i]]$grow_small <- yrBetas[[i]]$Intercept + yrBetas[[i]]$logarea.t0*size_info[[i]]$small ##-
#     #(size_info[[i]]$small*size_info[[i]]$scale + size_info[[i]]$center) # subtract rescaled focal size
#   yrBetas[[i]]$grow_big <- yrBetas[[i]]$Intercept + yrBetas[[i]]$logarea.t0*size_info[[i]]$large #-
#     #(size_info[[i]]$large*size_info[[i]]$scale + size_info[[i]]$center) # subtract rescaled focal size
# 
#   #merge in climate covariates
#   yrBetas[[i]] <- merge(yrBetas[[i]],Cdat,all.x=T)
# 
# }

### 1. fit INLA growth models with year and treatment effects ------------------------

# read in distance weights
dists <- read.csv(paste0(root,"/ExperimentTests/data/idaho_modern/speciesdata/IdahoModDistanceWeights_noExptl.csv"))

# import climate covariates
Cdat <- readRDS('data/temp_data/all_clim_covs.RDS')

# object to save year effects
sppList <- c("ARTR","HECO","POSE","PSSP")
yrBetas <- vector("list",length(sppList))
size_info <- vector("list",length(sppList))

for(i in 1:length(sppList)){

  iSpp <- sppList[i]

  # fit demography model
  source("analysis/growth/inla_growth.r")

  # get size info
  size_small <-quantile(allD$logarea.t0,0.1) ; size_large <- quantile(allD$logarea.t0,0.9) # get 10% and 90% size quantiles
  size_info[[i]] <- list( size_small, size_large)
  names(size_info[[i]]) <- c("small","large")
  rm(allD)

  #extract parameters for controls
  yrBetas[[i]] <- ranef(m1)$year
  yrBetas[[i]][,1] <-  yrBetas[[i]][,1] + fixef(m1)[1] ; yrBetas[[i]][,2] <-  yrBetas[[i]][,2] + fixef(m1)[2]  # add overall means
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

  # calculated relative growth rates of small and large plants on arithmetic scale
  yrBetas[[i]]$grow_small <- yrBetas[[i]]$Intercept + yrBetas[[i]]$logarea.t0*size_info[[i]]$small - size_info[[i]]$small
  yrBetas[[i]]$grow_big <- yrBetas[[i]]$Intercept + yrBetas[[i]]$logarea.t0*size_info[[i]]$large - size_info[[i]]$large

  #merge in climate covariates
  yrBetas[[i]] <- merge(yrBetas[[i]],Cdat,all.x=T)

}

# a <- yrBetas[[1]]$Intercept; b <- yrBetas[[1]]$logarea.t0
# # plot annual growth by size
# x <- seq(-1,7,0.01)
# plot(x,runif(length(x),-2,5),type="n",xlim=c(-1,7),ylim=c(-1,7),xlab="Size t_0",ylab="Size t_1")
# for(j in 1:25){
#   abline(a[j],b[j])
# }
# abline(0,1,col="red")

### 2. fit elastic net for growth of small then large plants ------------------------------------------------------------

enet_predictions <- vector("list",length(sppList))
best_coefs_small <- matrix(NA,length(sppList),ncol(Cdat)-4) # warning: hardwired "4" on climate columns
colnames(best_coefs_small) <- names(Cdat[,5:ncol(Cdat)])
best_coefs_big <- best_coefs_small

for(i in 1:length(sppList)){
  
  doSpp <- sppList[i]
  
  # prepare training data
  trainD <- yrBetas[[i]] %>% filter(year < 2011, year > 1926) # drop first year because of NAs in covariates
  y_small <- trainD$grow_small
  y_big <- trainD$grow_big
  X <- trainD[,9:NCOL(trainD)]
  X_mean <- colMeans(X); X_sd <- apply(X,MARGIN=2,FUN="sd")
  X <- scale(X, center = X_mean, scale = X_sd)
  
  # prepare out of sample data
  newD <- yrBetas[[i]] %>% filter(year >= 2011)
  y_new_small <- newD$grow_small
  y_new_big <- newD$grow_big
  X_new <- newD[,9:NCOL(newD)]
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
  saveRDS(best_coefs_small,"analysis/growth/enet_best_coefs_small.RDS")
  saveRDS(best_coefs_big,"analysis/growth/enet_best_coefs_big.RDS")

}

# clean up
rm(list= ls()[!(ls() %in% c('best_coefs_big','best_coefs_small','enet_predictions','sppList','yrBetas','size_info'))]) 


### 3. Draw figures -------------------------------------------  

# loop through species and plot observed and predicted intercepts and slopes

pdf("figures/enet_growth_yr_effects.pdf",height=4,width=8)
par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1))

for(i in 1:length(sppList)){
  
  # do small plants
  trts <- enet_predictions[[i]]$trts
  keep <- grep("small",names(enet_predictions[[i]]))
  fig_dat <- enet_predictions[[i]][keep]
  names(fig_dat) <- sub("_small","",names(fig_dat)) # trim names
  obs_pred_fig(fig_dat, trts, vital_rate="growth", effect="small plants", legend_location= "topleft")
  
  # do big plants
  trts <- enet_predictions[[i]]$trts
  keep <- grep("big",names(enet_predictions[[i]]))
  fig_dat <- enet_predictions[[i]][keep]
  names(fig_dat) <- sub("_big","",names(fig_dat)) # trim names
  obs_pred_fig(fig_dat, trts, vital_rate="growth", effect="big plants", legend_location= "topleft")
  
}

dev.off()















