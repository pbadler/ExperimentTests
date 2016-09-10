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

scale_clim_covs <- function(df , train, hold, clim_vars ){ 
  
  training_df <- df[train, ]
  holding_df  <- df[hold, ]
  
  clim_means <- colMeans(training_df[ , clim_vars], na.rm = TRUE)
  clim_sds   <- apply(training_df[ , clim_vars ], 2, FUN = sd, na.rm = TRUE)
  
  training_df[ , clim_vars] <- scale( training_df[ , clim_vars], center = TRUE, scale = TRUE)
  holding_df[  , clim_vars] <- scale( holding_df[  , clim_vars], center = clim_means, scale = clim_sds )
  
  return( list( training_df, holding_df) ) 
} 


make_datalist <- function(df, train, hold, clim_vars){ 
  
  # Function simply makes list of data for STAN models  
  
  # --------split into training and holding data and scale climate covariates -------------------------
  
  out <- scale_clim_covs(df, train, hold, clim_vars)
  
  training_df <- out[[1]]
  holding_df <- out[[2]]
  
  # --------training data -----------------------------------------------------------------------------
  N         <- nrow(training_df)                                  # number of data points for training data 
  nyrs      <- length(unique(training_df$year))                   # number of years 
  yid       <- as.numeric(as.factor(training_df$year))            # integer id for each year 
  Y         <- training_df$logarea.t1                             # plant size at time t 
  X         <- training_df$logarea.t0                             # plant size at time t - 1  
  
  C         <- as.matrix(training_df[, clim_vars])                # all climate covariates 
  Covs      <- ncol(C)                                            # number of climate covariates 
  
  W         <- as.matrix(training_df[, grep('W', names(training_df)) ]) # crowding matrix 
  W_covs    <- ncol(W)                                            # number of species in crowding matrix 
  
  gid       <- as.numeric(training_df$Group)                      # integer id for each plot area   
  G         <- length(unique(training_df$Group))                  # number of groups representing exclosure areas
  
  #---------hold out/prediction data ------------------------------------------------------------------
  npreds    <- nrow(holding_df)                                   # total predicted observations, modern data 
  y_holdout <- holding_df$logarea.t1                              # plant size at time t, modern data 
  Xhold     <- holding_df$logarea.t0                              # plant size at time t-1, modern data 
  Chold     <- holding_df[ , clim_vars]                           # climate matrix, modern data 
  Whold     <- holding_df[ , grep('W', names(holding_df)) ]       # crowding matrix, modern data 
  gid_out   <- as.numeric(holding_df$Group)                       # group id, modern data 
  
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

gf <- dir('data/temp_data/', pattern = 'growth.RDS', full.names = TRUE)
sf <- dir('data/temp_data/', pattern = 'survival.RDS', full.names = TRUE)
rf <- dir('data/temp_data/', pattern = 'recruitment.RDS', full.names = TRUE)

spp_names <- as.character( regmatches(c(gf, sf, rf), m = gregexpr( pattern = '([A-Z]{4})', c(gf, sf, rf) )) )

# -- read growth records -----------------------------------------------------------------#

growth <- lapply( gf, readRDS)

# -- merge with all survival records -----------------------------------------------------#


# -- merge with all recruitment records --------------------------------------------------#

growth <- lapply( growth, function(x){ subset(x, year > 1926 & !Treatment %in% c('No_shrub', 'No_grass'))} )

growth <- lapply(growth, merge, y = clim_covs, by = c('Treatment', 'Period', 'year')) 

# -- make training and holding subsets ----------------------------------------------------# 

training <- lapply( growth, function(x) { which(x$Period == "Historical" & x$Treatment == 'Control') } ) 
holding  <- lapply( growth, function(x) { which(x$Period == "Modern" )  } ) 

# -- prepare for stan ---------------------------------------------------------------------# 

datalists <- mapply( FUN = make_datalist, df = growth, train = training, hold = holding, MoreArgs = list( 'clim_vars' = clim_vars ), SIMPLIFY = FALSE)

names(datalists) <- spp_names

# ---- output ------------------------------------------------------------------------------# 

saveRDS( datalists, 'data/temp_data/data_lists_for_stan_models.RDS')


