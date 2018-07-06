
rm(list = ls() ) 

setwd("C:/Repos/ExperimentTests/precip/")
root <- "c:/repos/" # this is needed to get to the driversdata directory

#library(texreg) # to save output
library(xtable)
library(lme4)
library(dplyr)
library(glmnet) # package for statistical regularization
library(INLA)

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
  source("analysis/survival/inla_survival.r")
   
  #extract parameters for controls
  Intercept=m1$summary.random$yearID$mean
  year=m1$summary.random$year$ID
  logarea=m1$summary.random$year$mean
  yrBetas[[i]]=data.frame(year,Intercept,logarea)
  yrBetas[[i]]$Treatment <- "Control"
  yrBetas[[i]]$year <-  as.numeric(as.character(yrBetas[[i]]$year))

  #add parameters for drought treatment
  tmp <- subset(yrBetas[[i]],year > 2010 & Treatment == "Control")
  tmp$Treatment <- "Drought"
  tmp$Intercept <- tmp$Intercept + m1$summary.fixed$mean[which(m1$names.fixed=="TreatmentDrought")]
  yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
  
  #add parameters for irrigation treatment
  tmp <- subset(yrBetas[[i]], year > 2010 & Treatment == "Control")
  tmp$Treatment <- "Irrigation"
  tmp$Intercept <- tmp$Intercept + m1$summary.fixed$mean[which(m1$names.fixed=="TreatmentIrrigation")]
  yrBetas[[i]] <- rbind(yrBetas[[i]] ,tmp)
  
  #merge in climate covariates
  yrBetas[[i]] <- merge(yrBetas[[i]],Cdat,all.x=T)
  
}

### elastic net ------------------------------------------------------------

best_coefs <- matrix(NA,length(sppList),ncol(Cdat)-4) # warning: hardwired "4" on climate columns
colnames(best_coefs) <- names(Cdat[,5:ncol(Cdat)])
rmse_ratio <- numeric(length(sppList))

pdf("figures/enet_survival_yr_effects.pdf",height=4,width=4)
par(tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,2,1))

for(i in 1:length(sppList)){
  
  doSpp <- sppList[i]
  
  ## model Intercept year effects
  # prepare data
  trainD <- yrBetas[[i]] %>% filter(year < 2011, year > 1926) # drop first year because of NAs in covariates
  y <- trainD$Intercept
  X <- trainD[,7:NCOL(trainD)]
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
  
  # get out of sample MSE
  newD <- yrBetas[[i]] %>% filter(year >= 2011)
  y_new <- newD$Intercept
  X_new <- newD[,7:NCOL(newD)]
  X_new <- scale(X_new, center = X_mean, scale = X_sd)
  y_hat_new <- predict(enet_out,newx=X_new, s="lambda.min")
  mse_new <- mean((y_new-y_hat_new)^2)
  
  # make figure
  plot(c(y,y_new),c(y_hat,y_hat_new),type="n",xlab="Observed",ylab="Predicted",
       ylim=c(min(c(y,y_new,y_hat,y_hat_new)),max(c(y,y_new,y_hat,y_hat_new))),
       xlim=c(min(c(y,y_new,y_hat,y_hat_new)),max(c(y,y_new,y_hat,y_hat_new))),
       main=paste0(sppList[i]," survival year effects (Intercept)"))
  abline(0,1)
  points(y,y_hat)
  points(y_new[which(newD$Treatment=="Control")],y_hat_new[which(newD$Treatment=="Control")],pch=16)
  points(y_new[which(newD$Treatment=="Drought")],y_hat_new[which(newD$Treatment=="Drought")],pch=16,col="red")
  points(y_new[which(newD$Treatment=="Irrigation")],y_hat_new[which(newD$Treatment=="Irrigation")],pch=16,col="blue")
  legend("bottomright",c("Control (training)","Control (out-of-sample)","Drought (out-of-sample)",
                     "Irrigation (out-of-sample)"),pch=c(1,16,16,16),
                      col=c("black","black","red","blue"),bty="n",cex=0.8)
}

dev.off()

saveRDS(best_coefs,"analysis/survival/enet_best_coefs.RDS")
