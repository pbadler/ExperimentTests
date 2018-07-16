
rm(list = ls() ) 

setwd("C:/Repos/ExperimentTests/precip/")
root <- "c:/repos/" # this is needed to get to the driversdata directory

#library(texreg) # to save output
library(xtable)
library(lme4)
library(dplyr)
library(glmnet) # package for statistical regularization


### 1. import year effects from stan model ------------------------

# import climate covariates
Cdat <- readRDS('data/temp_data/all_clim_covs.RDS')

# object to save year effects
sppList <- c("ARTR","HECO","POSE","PSSP")
yrBetas <- vector("list",length(sppList))
size_info <- vector("list",length(sppList))

for(i in 1:length(sppList)){
 
  iSpp <- sppList[i]
  
  # get vector of years
  dat <- readRDS("data/temp_data/growth_data_lists_for_stan.RDS")[[i]]
  year <- unique(dat$year2)
  size_small <-quantile(dat$X,0.1) ; size_large <- quantile(dat$X,0.9) # get 10% and 90% size quantiles
  size_info[[i]] <- c( size_small, size_large, dat$Xcenter, dat$Xscale )
  names(size_info[[i]]) <- c("small","large","center","scale")
  rm(dat)
  
  # import model summary
  m1 <- read.csv(paste0("output/treatment_model_parameters_",iSpp,"_growth.csv"))
   
  #extract parameters for controls
  keep <- grep("a", substring(as.character(m1$X),1,1))
  Intercept <- m1$mean[keep]
  keep <- grep("b1", substring(as.character(m1$X),1,2))
  keep <- keep[m1$X[keep]!="b1_mu"] 
  logarea.t0 <- m1$mean[keep]
  yrBetas[[i]] <- data.frame(year,Intercept,logarea.t0)
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
  
  #merge in climate covariates
  yrBetas[[i]] <- merge(yrBetas[[i]],Cdat,all.x=T)
  
}

### 2. fit elastic net on intercept and slope ------------------------------------------------------------

enet_predictions <- vector("list",length(sppList))
best_coefs <- matrix(NA,length(sppList),ncol(Cdat)-4) # warning: hardwired "4" on climate columns
colnames(best_coefs) <- names(Cdat[,5:ncol(Cdat)])
best_coefs_slope <- best_coefs

for(i in 1:length(sppList)){
  
  doSpp <- sppList[i]
  
  # prepare training data
  trainD <- yrBetas[[i]] %>% filter(year < 2011, year > 1926) # drop first year because of NAs in covariates
  y_int <- trainD$Intercept
  y_slope <- trainD$logarea.t0
  X <- trainD[,7:NCOL(trainD)]
  X_mean <- colMeans(X); X_sd <- apply(X,MARGIN=2,FUN="sd")
  X <- scale(X, center = X_mean, scale = X_sd)
  
  # prepare out of sample data
  newD <- yrBetas[[i]] %>% filter(year >= 2011)
  y_new_int <- newD$Intercept
  y_new_slope <- newD$logarea.t0
  X_new <- newD[,7:NCOL(newD)]
  X_new <- scale(X_new, center = X_mean, scale = X_sd)
  
  # fit elastic net on intercept
  pen_facts <- rep(1, ncol(X)) # penalize all covariates
  lambdas <- 10^seq(2, -2, by = -0.005) # sequence of penalties to test
  enet_int <- cv.glmnet(x = X, 
                       y = y_int, 
                       lambda = lambdas,
                       penalty.factor = pen_facts,
                       family = "gaussian", 
                       alpha = 0.5, # 0 for ridge, 1 for lasso 
                       standardize = FALSE, 
                       type.measure = "mse",
                       nfolds = length(y_int))
  
  best_coefs[i,] <- enet_int$glmnet.fit$beta[,which(enet_int$lambda==enet_int$lambda.min)]
  
  # predictions for training data
  y_hat_int <- predict(enet_int,newx=X,s="lambda.min")
  
  # predictions for out of sample data
  y_hat_new_int <- predict(enet_int,newx=X_new, s="lambda.min")

  # fit elastic net on slope year effects
  enet_slope <- cv.glmnet(x = X, 
                       y = y_slope, 
                       lambda = lambdas,
                       penalty.factor = pen_facts,
                       family = "gaussian", 
                       alpha = 0.5, # 0 for ridge, 1 for lasso 
                       standardize = FALSE, 
                       type.measure = "mse",
                       nfolds = length(y_slope))
  
  best_coefs_slope[i,] <- enet_slope$glmnet.fit$beta[,which(enet_slope$lambda==enet_slope$lambda.min)]
  
  # predictions for training data
  y_hat_slope <- predict(enet_slope,newx=X,s="lambda.min")
  
  # predictions for out of sample data
  y_hat_new_slope <- predict(enet_slope,newx=X_new, s="lambda.min")
  
  # save output
  enet_predictions[[i]] <- list(trts = newD$Treatment,
                                y_int=y_int,y_hat_int=y_hat_int,y_new_int=y_new_int, y_hat_new_int=y_hat_new_int,
                                y_slope=y_slope, y_hat_slope=y_hat_slope, y_new_slope=y_new_slope, y_hat_new_slope=y_hat_new_slope)
  saveRDS(best_coefs,"analysis/growth/enet_best_coefs_intercept.RDS")
  saveRDS(best_coefs_slope,"analysis/growth/enet_best_coefs_slope.RDS")

}

# clean up
rm(list= ls()[!(ls() %in% c('best_coefs','best_coefs_slope','enet_predictions','sppList','yrBetas','size_info'))]) 


### 3. Draw figures -------------------------------------------  

# custom plotting function
obs_pred_fig <- function(my_dat, vital_rate, effect="int", legend_location){
  # my_dat is a list containing observations and predictions
  # vital_rate is a string: either "growth" or "survival"
  # effect is a string with default "int" and alternative value "slope"
  # legend_location is a string, such as "topleft"
  
  # subset data of interest
  trts <- my_dat$trts
  keep <- grep(effect,names(my_dat))
  my_dat <- my_dat[keep]
  
  # trim names
  names(my_dat) <- sub(paste0("_",effect),"",names(my_dat))
  
  attach(my_dat)
  plot(c(y,y_new),c(y_hat,y_hat_new),type="n",xlab="Observed",ylab="Predicted",
       ylim=c(min(c(y,y_new,y_hat,y_hat_new)),max(c(y,y_new,y_hat,y_hat_new))),
       xlim=c(min(c(y,y_new,y_hat,y_hat_new)),max(c(y,y_new,y_hat,y_hat_new))),
       main=paste(sppList[i],vital_rate,effect, "year effects"))
  abline(0,1)
  points(y,y_hat)
  points(y_new[which(trts =="Control")],y_hat_new[which(trts =="Control")],pch=16)
  points(y_new[which(trts =="Drought")],y_hat_new[which(trts =="Drought")],pch=16,col="red")
  points(y_new[which(trts =="Irrigation")],y_hat_new[which(trts =="Irrigation")],pch=16,col="blue")
  legend(legend_location,c("Control (training)","Control (out-of-sample)","Drought (out-of-sample)",
                     "Irrigation (out-of-sample)"),pch=c(1,16,16,16),
                      col=c("black","black","red","blue"),bty="n",cex=0.8)
  detach(my_dat)
}

# loop through species and plot observed and predicted intercepts and slopes

pdf("figures/enet_growth_yr_effects.pdf",height=4,width=8)
par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1))

for(i in 1:length(sppList)){
  
  obs_pred_fig(enet_predictions[[i]], vital_rate="growth", effect="int", legend_location= "topleft")
  obs_pred_fig(enet_predictions[[i]], vital_rate="growth", effect="slope", legend_location= "topleft")
  
}

dev.off()














# ### 1. fit growth models with year and treatment effects ------------------------
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
#   source("analysis/growth/inla_growth.r")
#    
#   #extract parameters for controls
#   yrBetas[[i]] <- ranef(m1)$year
#   yrBetas[[i]]$year <- row.names(yrBetas[[i]])
#   yrBetas[[i]]$Treatment <- "Control"
#   names(yrBetas[[i]])[1] <- "Intercept"
#   
#   #add parameters for drought treatment
#   tmp <- subset(yrBetas[[i]],year > 2010)
#   tmp$Treatment <- "Drought"
#   tmp$Intercept <- tmp$Intercept + fixef(m1)[which(names(fixef(m1))=="TreatmentDrought")]
#   yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
#   
#   #add parameters for irrigation treatment
#   tmp <- subset(yrBetas[[i]],year > 2010 & Treatment=="Control")
#   tmp$Treatment <- "Irrigation"
#   tmp$Intercept <- tmp$Intercept + fixef(m1)[which(names(fixef(m1))=="TreatmentIrrigation")]
#   yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
#   
#   #merge in climate covariates
#   yrBetas[[i]] <- merge(yrBetas[[i]],Cdat,all.x=T)
#   
# }



