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

scale_clim_covs <- function( df, train, hold, clim_vars ){ 
  
  clim_means <- colMeans(df[ train, clim_vars], na.rm = TRUE)
  clim_sds   <- apply(df[ train, clim_vars ], 2, FUN = sd, na.rm = TRUE)
  
  df[train, clim_vars] <- scale( df[ train, clim_vars], center = TRUE, scale = TRUE)
  df[hold, clim_vars]  <- scale( df[ hold, clim_vars ], center = clim_means, scale = clim_sds )
  
  df
} 


make_datalist <- function(df, train, hold, clim_vars){ 
  
  # Function simply makes list of data for STAN models  
  
  # --------split into training and holding data ------------------------------------------------------
  
  training_df <- df[train, ]
  holding_df  <- df[hold, ]
  
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

growth_files <- dir('data/temp_data/', pattern = 'growth.RDS', full.names = TRUE)

#survival_files <- dir('data/temp_data/')
# recruitment_files <- dir('data/temp_data/)

spp_names <- regmatches(growth_files, m = gregexpr( pattern = '([A-Z]{4})', growth_files )) 

# -- merge with all growth records -------------------------------------------------------#

growth <- lapply( growth_files, readRDS)

growth <- lapply( growth, function(x){ subset(x, year > 1926 & !Treatment %in% c('No_shrub', 'No_grass'))} )

growth <-  lapply(growth, merge, y = clim_covs, by = c('Treatment', 'Period', 'year')) 

# -- merge with all survival records -----------------------------------------------------#


# -- merge with all recruitment records --------------------------------------------------#


# -- scale climate based on historical period ------------------------------------------# 

training <- lapply( growth, function(x) { which(x$Period == "Historical" & x$Treatment == 'Control') } ) 
holding  <- lapply( growth, function(x) { which(x$Period == "Modern" )  } ) 

growth <- mapply( FUN = scale_clim_covs, df = growth, train = training, hold = holding, MoreArgs = list( 'clim_vars' = clim_vars )) 
# ggplot(subset( test[, c('year', 'Treatment', 'Period', clim_vars)], Period == "Historical") %>% gather_('var', 'val' , clim_vars) , aes( x = year, y = val, color = Treatment )) + geom_point() + geom_line() +  facet_wrap(~ var ) 
# ggplot(subset( test[, c('year', 'Treatment', 'Period', clim_vars)], Period == "Modern") %>% gather_('var', 'val' , clim_vars) , aes( x = year, y = val, color = Treatment )) + geom_point() + geom_line() +  facet_wrap(~ var ) 

# prepare for stan -----------------------------------------------------------------------# 

datalist <- mapply( FUN = make_datalist, df = growth, train = training, hold = holding, MoreArgs = list( 'clim_vars' = clim_vars ), SIMPLIFY = FALSE)

datalist[[1]]$yid

datalist[[1]]$yid

# ---- output ----------------------------------------------------------------------------# 

saveRDS( datalist, 'data/temp_data/data_lists_for_stan_models.RDS')


