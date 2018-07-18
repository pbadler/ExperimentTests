
rm(list = ls() ) 

setwd("C:/Repos/ExperimentTests/precip/")
root <- "c:/repos/" # this is needed to get to the driversdata directory

#library(texreg) # to save output
library(xtable)
library(lme4)
library(dplyr)
library(glmnet) # package for statistical regularization
library(INLA)
source("analysis/figure_scripts/elastic_net_observed_vs_predicted.R")

### 1. Import year effects from Stan models -------------------------------

# import climate covariates
Cdat <- readRDS('data/temp_data/all_clim_covs.RDS')

# object to save year effects
sppList <- c("ARTR","HECO","POSE","PSSP")
yrBetas <- vector("list",length(sppList))

for(i in 1:length(sppList)){
 
  iSpp <- sppList[i]
  
  # get vector of years
  dat <- readRDS("data/temp_data/recruitment_data_lists_for_stan.RDS")[[i]]
  year <- unique(dat$year2)
  rm(dat)
  
  # import model summary
  m1 <- read.csv(paste0("output/treatment_model_parameters_",iSpp,"_recruitment.csv"))
   
  #extract parameters for controls

  keep <- grep("a", substring(as.character(m1$X),1,1))
  Intercept=m1$mean[keep]
  yrBetas[[i]]=data.frame(year,Intercept)
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

### 2. elastic net on intercept year effects ------------------------------------------------------------

best_coefs <- matrix(NA,length(sppList),ncol(Cdat)-4) # warning: hardwired "4" on climate columns
colnames(best_coefs) <- names(Cdat[,5:ncol(Cdat)])
rmse_ratio <- numeric(length(sppList))

pdf("figures/enet_recruit_yr_effects.pdf",height=4,width=4)
par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1))

for(i in 1:length(sppList)){
  
  doSpp <- sppList[i]
  
  ## model Intercept year effects
  # prepare data
  trainD <- yrBetas[[i]] %>% filter(year < 2011)
  y <- trainD$Intercept
  X <- trainD[,6:NCOL(trainD)]
  X_mean <- colMeans(X); X_sd <- apply(X,MARGIN=2,FUN="sd")
  X <- scale(X, center = X_mean, scale = X_sd)
  pen_facts <- rep(1, ncol(X)) # penalize all covariates
  lambdas <- 10^seq(2, -2, by = -0.005) # sequence of penalties to test
  
  
  enet_out <- cv.glmnet(x = X, 
                       y = y, 
                       lambda = lambdas,
                       penalty.factor = pen_facts,
                       family = "gaussian", 
                       alpha = 0.5, # 0 for ridge, 1 for lasso 
                       standardize = FALSE, 
                       type.measure = "mse",
                       nfolds = length(y))
  
  # look at results
  # plot(log(enet_out$lambda),enet_out$cvm,xlab="log(Lambda)",ylab="CV score",type="l")
  
  #matplot(log(enet_out$lambda),t(enet_out$glmnet.fit$beta),type="l")
  #abline(v=log(enet_out$lambda.min),lty="dashed",col="darkgrey")
  
  best_coefs[i,] <- enet_out$glmnet.fit$beta[,which(enet_out$lambda==enet_out$lambda.min)]
  
  # predictions for training data
  y_hat <- predict(enet_out,newx=X,s="lambda.min")
  
  # get out of sample predictions
  newD <- yrBetas[[i]] %>% filter(year >= 2011)
  y_new <- newD$Intercept
  X_new <- newD[,6:NCOL(newD)]
  X_new <- scale(X_new, center = X_mean, scale = X_sd)
  y_hat_new <- predict(enet_out,newx=X_new, s="lambda.min")
  trts <- newD$Treatment
  
  # exponentiate to get per capita recruitment rate
  y <-exp(y); y_hat <- exp(y_hat)
  y_new <- exp(y_new); y_hat_new <- exp(y_hat_new)
  
  
  fig_dat <- list(y=y, y_hat=y_hat, y_new=y_new, y_hat_new=y_hat_new)
  obs_pred_fig(fig_dat, trts, vital_rate="recruitment", effect="", legend_location= "topleft")
  
}

dev.off()

saveRDS(best_coefs,"analysis/survival/enet_best_coefs_intercept.RDS")

