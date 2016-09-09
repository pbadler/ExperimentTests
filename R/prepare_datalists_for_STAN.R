######################################################################################
#
# Make STAN datalist  
#
#####################################################################################

rm(list = ls() )

library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)


make_datalist <- function(train, hold){ 
  
  # Function simply makes list of data for STAN models  

  # --------training data -----------------------------------------------------------------------------
  N         <- nrow(train)                                  # number of data points for training data 
  nyrs      <- length(unique(train$year))                   # number of years 
  yid       <- as.numeric(as.factor(train$year))            # integer id for each year 
  Y         <- train$logarea.t1                             # plant size at time t 
  X         <- train$logarea.t0                             # plant size at time t - 1  
  
  C         <- as.matrix(train[, clim_vars])                # all climate covariates 
  Covs      <- ncol(C)                                      # number of climate covariates 
  
  W         <- as.matrix(train[, grep('W', names(train)) ]) # crowding matrix 
  W_covs    <- ncol(W)                                      # number of species in crowding matrix 
  
  gid       <- as.numeric(train$Group)                      # integer id for each plot area   
  G         <- length(unique(train$Group))                  # number of groups representing exclosure areas
  
  #---------hold out/prediction data ------------------------------------------------------------------
  npreds    <- nrow(hold)                                   # total predicted observations, modern data 
  y_holdout <- hold$logarea.t1                              # plant size at time t, modern data 
  Xhold     <- hold$logarea.t0                              # plant size at time t-1, modern data 
  Chold     <- hold[ , clim_vars]                           # climate matrix, modern data 
  Whold     <- hold[ , grep('W', names(hold)) ]             # crowding matrix, modern data 
  gid_out   <- as.numeric(hold$Group)                       # group id, modern data 
  
  return( 
    list(
      N = N, Y = Y, X = X , gid = gid, G = G, Yrs = nyrs, yid = yid , Covs = Covs, C = C, W = W, W_covs = W_covs,
      gid_out = gid_out, npreds = npreds, y_holdout = y_holdout, Xhold = Xhold, Chold = Chold, Whold = Whold, 
      tau_beta = 1
    )
  )
}



# -- select covariates -------------------------------------------------------------------#

clim_vars <- c('T.sp.1', 'T.sp.2', 'P.w.sp.1', 'P.w.sp.2', 'T.su.1', 'T.su.2', 'P.a.0', 'P.a.1')
clim_vars <- sort(clim_vars)

# -- read data files ---------------------------------------------------------------------# 

clim_covs <- readRDS('data/temp_data/all_clim_covs.RDS')

growth_files <- dir('data/temp_data/', pattern = 'growth.RDS', full.names = TRUE)

#survival_files <- dir('data/temp_data/')
# recruitment_files <- dir('data/temp_data/)

spp_names <- regmatches(growth_files, m = gregexpr( pattern = '([A-Z]{4})', growth_files )) 

# -- merge with all growth records -------------------------------------------------------#

growth <- lapply( growth_files, readRDS)

growth <-  lapply(growth, left_join, y = clim_covs, by = c('Treatment', 'Period', 'year')) 

names(growth) <- spp_names

# -- merge with all survival records -----------------------------------------------------#


# -- merge with all recruitment records --------------------------------------------------#


# prepare for stan -----------------------------------------------------------------------# 
test <- growth[[1]] # test one species 

# use historical data as training data ------------------------------------------------# 
train <- subset(test, Period == 'Historical' & year > 1926) 

# use modern data as hold out data ----------------------------------------------------#
hold <- subset(test, Period == 'Modern' & ! Treatment %in% c('No_grass', 'No_shrub') )

datalist <- make_datalist(train, hold)

# ---- output ----------------------------------------------------------------------------# 

saveRDS( datalist, 'data/temp_data/data_lists_for_stan_models.RDS')


