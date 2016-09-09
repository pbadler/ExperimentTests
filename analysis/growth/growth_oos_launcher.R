##
##  R script for growth as run on Utah State University's
##  High Performance Computing system.
## 
##  Adapted from Andrew Tredennick's script 
##
##  Script takes command line argument 'do_grid' to run all levels of
##  regularization and each cross-validation set in parallel. Script launches
##  n*k models: 'n' prior standard deviations and 'k' validation sets.
##
##  Saved output:
##    1. Summed log pointwise predictive density
##    2. WAIC score which approximates out of sample predictive score
##
##  Original Author:      Andrew Tredennick
##  Modified By:          Andrew Kleinhesselink 
##  Email:                arklein@aggiemail.usu.edu
##  Date created:         12-6-2015
##  Date modified:        09-08-2016
##

# Change species four letter code here...

# -- read in species data --------------------------------------------------------# 

G <- readRDS('data/temp_data/all_growth_combined.RDS')

test <- G[[1]] # test one species 

test <- subset(test, ! Treatment %in% c('No_grass', 'No_shrub'))

test <- test[order(test$year), ]

####
####  Set SD Prior and CV Set from Command Line Arguments
####

args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_grid <- as.numeric(myargument)

####
####  Load Libraries and Subroutines
####

library(rstan)
# library(plyr)
# library(reshape2)
library(ggmcmc)
library(matrixStats)
source("waic_fxns.R")

####
####  Compile Stan Model
####

#- scale covariates ----------------------------------------------------------------------------# 

clim_vars <- c('T.sp.1', 'T.sp.2', 'P.w.sp.1', 'P.w.sp.2', 'T.su.1', 'T.su.2', 'P.a.0', 'P.a.1')
clim_vars <- sort(clim_vars)

# make interactions: 
# test$TxP.sp.1 <- test$P.sp.1*test$T.sp.1
# test$TxP.sp.2 <- test$P.sp.2*test$T.sp.2

#df_list <- lapply( df_list, function( x ) { x$ppt})

clim_covs_1 <- test[ test$period == 'historical'  , clim_vars ]
clim_covs_2 <- test[ test$period == 'contemporary', clim_vars ]

# Get scalers for climate covariates from historical data
clim_means <- colMeans(clim_covs_1, na.rm = T)
clim_sds <- apply(clim_covs_1, 2, FUN = sd, na.rm = TRUE)

clim_covs_1 <- scale(clim_covs_1, center = TRUE, scale = TRUE) # scale historical data 

# scale contemporary and treatment data by historical mean and sd 

head(clim_covs_1)
head(clim_covs_2)
tail(clim_covs_2)
clim_covs_2 <- 
head( scale( clim_covs_2, center = clim_means, scale = clim_sds ) )

head( clim_covs ) 

groups <- as.numeric(trainD$Group)
G <- length(unique(trainD$Group))
nyrs <- length(unique(trainD$year))
W <- cbind(trainD$W, trainD$W*log(trainD$area.t0))
yid <- as.numeric(as.factor(trainD$year))

# Holdout data
holdD <- subset(growD, year==holdyear)
holdD$ppt1TmeanSpr1 <- holdD$ppt1*holdD$TmeanSpr1
holdD$ppt2TmeanSpr2 <- holdD$ppt2*holdD$TmeanSpr2
# holdD$sizepptLag <- holdD$pptLag*log(holdD$area.t0)
# holdD$sizeppt1 <- holdD$ppt1*log(holdD$area.t0)
# holdD$sizeppt2 <- holdD$ppt2*log(holdD$area.t0)
# holdD$sizeTmeanSpr1 <- holdD$TmeanSpr1*log(holdD$area.t0)
# holdD$sizeTmeanSpr2 <- holdD$TmeanSpr2*log(holdD$area.t0)
clim_covs_oos <- holdD[,clim_vars_all]
for(j in 1:ncol(clim_covs_oos)){
  clim_covs_oos[,j] <- (clim_covs_oos[,j] - clim_means[j])/clim_sds[j]
}
W_oos <- cbind(holdD$W, holdD$W*log(holdD$area.t0))
gid_out <- as.numeric(holdD$Group) 

datalist <- list(N=nrow(trainD), Yrs=nyrs, yid=yid,
                 Covs=ncol(clim_covs), Y=log(trainD$area.t1), X=log(trainD$area.t0),
                 C=clim_covs, W=W, G=G, gid=groups, tau_beta=1,
                 npreds=nrow(holdD), y_holdout=log(holdD$area.t1), Xhold=log(holdD$area.t0),
                 Chold=clim_covs_oos, Whold=W_oos, gid_out=gid_out)
pars <- c("log_lik")
mcmc_oos <- stan(file="growth_oos_cv.stan", data=datalist, 
                 pars=pars, chains=0)


####
####  Set up CV x Regularization grid
####
n.beta <- 24
sd_vec <- seq(0.1,1.5,length.out = n.beta)
yrs.vec <- unique(growD$year)
K <- length(yrs.vec)
cv.s2.grid <- expand.grid(1:n.beta,1:K)
n.grid <- dim(cv.s2.grid)[1]
fold.idx.mat <- matrix(1:length(yrs.vec),ncol=K)



####
####  Source Growth Model Function and Fit Model
####
source("growth_fxns.R")
out_lpd <- cv.fcn(do_grid)
saveRDS(out_lpd, paste0(do_species,"_oos_cv_dogrid_",do_grid,".RDS"))


