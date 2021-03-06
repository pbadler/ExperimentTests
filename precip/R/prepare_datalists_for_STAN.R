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
library(stringr)

scale_clim_covs <- function(df , train, hold, clim_vars ){ 
  
  training_df <- df[train, ]
  holding_df  <- df[hold, ]
  
  clim_means <- colMeans(training_df[ , clim_vars], na.rm = TRUE)
  clim_sds   <- apply(training_df[ , clim_vars ], 2, FUN = sd, na.rm = TRUE)
  
  training_df[ , clim_vars] <- scale( training_df[ , clim_vars], center = TRUE, scale = TRUE)
  holding_df[  , clim_vars] <- scale( holding_df[  , clim_vars], center = clim_means, scale = clim_sds )
  
  training_df <-  training_df [ order(training_df$year), ]
  holding_df  <-  holding_df [ order(holding_df$year), ]
  
  return( list( training_df, holding_df) ) 
} 


growth_dataframe2datalist <- function(df, train, hold, clim_vars){ 
  
  # Function simply makes list of data for STAN models  
  
  survival_dlist <- readRDS('data/temp_data/survival_data_lists_for_stan.RDS')  
    
  # --------split into training and holding data and scale climate covariates -------------------------
  
  out <- scale_clim_covs(df, train, hold, clim_vars)
  
  training_df <- out[[1]]
  holding_df <- out[[2]]
  
  # --------training data -----------------------------------------------------------------------------
  N         <- nrow(training_df)                                  # number of data points for training data 
  nyrs      <- length(unique(training_df$year))                   # number of years 
  yid       <- as.numeric(as.factor(training_df$year))            # integer id for each year 
  Y         <- training_df$logarea.t1                             # plant size at time t + 1 
  X         <- training_df$logarea.t0                             # plant size at time t   
  
  C         <- as.matrix(training_df[, clim_vars])                # all climate covariates 
  Covs      <- ncol(C)                                            # number of climate covariates 
  
  W         <- as.matrix(training_df[, grep('W', names(training_df)) ]) # crowding matrix 
  W_covs    <- ncol(W)                                            # number of species in crowding matrix 
  
  gid       <- as.numeric(training_df$Group)                      # integer id for each plot area   
  G         <- length(unique(training_df$Group))                  # number of groups representing exclosure areas
  
  trackid   <- training_df$trackID
  year      <- training_df$year
  quad      <- training_df$quad
  
  #---------hold out/prediction data ------------------------------------------------------------------
  npreds    <- nrow(holding_df)                                   # total predicted observations, modern data 
  y_holdout <- holding_df$logarea.t1                              # plant size at time t, modern data 
  Xhold     <- holding_df$logarea.t0                              # plant size at time t-1, modern data 
  Chold     <- holding_df[ , clim_vars]                           # climate matrix, modern data 
  Whold     <- holding_df[ , grep('W', names(holding_df)) ]       # crowding matrix, modern data 
  gid_out   <- as.numeric(holding_df$Group)                       # group id, modern data 
  yid_out   <- as.numeric(as.factor(holding_df$year))             # year id, modern data
  nyrs_out  <- length(unique(holding_df$year))                    # num years, modern data 
  treat_out <- as.numeric(factor(holding_df$Treatment))           # treatments 
  
  trackid_out   <- holding_df$trackID
  year_out      <- holding_df$year
  quad_out      <- holding_df$quad
  
  #--------all survival data for cover predictions ------------------------------------------------------------# 
  species_name   <- unique(training_df$species )
  spp_list  <- factor( str_extract( colnames(W), '[A-Z]+$'))       # get species names 
  spp       <- grep(species_name, spp_list)                       # assign species number to datalist 
  
  survival_df <- survival_dlist[[spp]]                             # get frame for the right species
  
  N2        <- nrow(survival_df)                                   # total predictions
  X2        <- survival_df$logarea.t0                              # plant size at time t-1 
  C2        <- survival_df$C                                        # climate matrix 
  W2        <- survival_df$W                                        # crowding matrix 
  gid2      <- as.numeric(survival_df$Group)                       # group id 
  yid2      <- as.numeric(as.factor(survival_df$year))             # year id
  nyrs2     <- length(unique(survival_df$year))                    # num years 
  treat2    <- as.numeric(factor(survival_df$Treatment))           # treatments 
  
  trackid2   <- survival_df$trackID
  year2      <- survival_df$year
  quad2      <- survival_df$quad
  
  return( 
    list(
      N = N, Y = Y, X = X , gid = gid, G = G, Yrs = nyrs, yid = yid , Covs = Covs, C = C, W = W, W_covs = W_covs,
      gid_out = gid_out, npreds = npreds, y_holdout = y_holdout, Xhold = Xhold, Chold = Chold, Whold = Whold, yid_out = yid_out, nyrs_out = nyrs_out, treat_out = treat_out,
      trackid = trackid, trackid_out = trackid_out, year = year, year_out = year_out, quad = quad, quad_out = quad_out,        
      tau_beta = 1, 
      spp = spp, 
      N2 = N2, X2 = X2, C2 = C2, W2 = W2, gid2 = gid2, yid2 = yid2, nrys2 = nyrs2, treat2 = treat2, trackid2 = trackid2, 
      year2 = year2, quad2 = quad2 
    )
  )
}

survival_dataframe2datalist <- function(df, train, hold, clim_vars){

  # Function simply makes list of data for STAN models

  # --------split into training and holding data and scale climate covariates -------------------------

  out <- scale_clim_covs(df, train, hold, clim_vars)

  training_df <- out[[1]]
  holding_df <- out[[2]]

  # --------training data -----------------------------------------------------------------------------
  N         <- nrow(training_df)                                  # number of data points for training data
  nyrs      <- length(unique(training_df$year))                   # number of years
  yid       <- as.numeric(as.factor(training_df$year))            # integer id for each year
  Y         <- training_df$survives                               # plant survival at time t + 1   
  X         <- training_df$logarea                                # plant size at time t

  C         <- as.matrix(training_df[, clim_vars])                # all climate covariates
  Covs      <- ncol(C)                                            # number of climate covariates

  W         <- as.matrix(training_df[, grep('W', names(training_df)) ]) # crowding matrix
  W_covs    <- ncol(W)                                            # number of species in crowding matrix

  gid       <- as.numeric(training_df$Group)                      # integer id for each plot area
  G         <- length(unique(training_df$Group))                  # number of groups representing exclosure areas

  trackid   <- training_df$trackID
  year      <- training_df$year
  quad      <- training_df$quad
  
  #---------hold out/prediction data ------------------------------------------------------------------
  species_name   <- unique(training_df$species )
  spp_list  <- factor( str_extract( colnames(W), '[A-Z]+$'))      # get species names 
  spp       <- grep(species_name, spp_list)                       # assign species number to datalist 
  
  npreds    <- nrow(holding_df)                                   # total predicted observations, modern data
  y_holdout <- holding_df$survives                                # plant survival at time t + 1, modern data
  Xhold     <- holding_df$logarea                                 # plant size at time t, modern data
  Chold     <- holding_df[ , clim_vars]                           # climate matrix, modern data
  Whold     <- holding_df[ , grep('W', names(holding_df)) ]       # crowding matrix, modern data
  gid_out   <- as.numeric(holding_df$Group)                       # group id, modern data
  yid_out   <- as.numeric(as.factor(holding_df$year))             # year id, modern data
  nyrs_out  <- length(unique(holding_df$year))                    # num years, modern data
  treat_out <- as.numeric(factor(holding_df$Treatment))           # treatments

  trackid_out   <- holding_df$trackID
  year_out      <- holding_df$year
  quad_out      <- holding_df$quad
  
  return(
    list(
      N = N, Y = Y, X = X , gid = gid, G = G, Yrs = nyrs, yid = yid , Covs = Covs, C = C, W = W, W_covs = W_covs,
      gid_out = gid_out, npreds = npreds, y_holdout = y_holdout, Xhold = Xhold, Chold = Chold, Whold = Whold, yid_out = yid_out, nyrs_out = nyrs_out, treat_out = treat_out,
      trackid = trackid, trackid_out = trackid_out, year = year, year_out = year_out, quad = quad, quad_out = quad_out,              
      tau_beta = 1, 
      spp = spp
    )
  )
}



recruitment_dataframe2datalist <- function(df, train, hold, clim_vars){
  
  # Function simply makes list of data for STAN models
  
  # --------split into training and holding data and scale climate covariates -------------------------
  
  out <- scale_clim_covs(df, train, hold, clim_vars)
  
  training_df <- out[[1]]
  holding_df <- out[[2]]
  
  # --------training data -----------------------------------------------------------------------------
  N         <- nrow(training_df)                                  # number of data points for training data
  
  nyrs      <- length(unique(training_df$year))                   # number of years
  yid       <- as.numeric(as.factor(training_df$year))            # integer id for each year
  Y         <- training_df$Y                                      # new recruit at time t + 1   
  
  C         <- as.matrix(training_df[, clim_vars])                # all climate covariates
  Covs      <- ncol(C)                                            # number of climate covariates
  
  parents1  <- as.matrix(training_df[, grep('^cov', names(training_df)) ])/100  # parents in plot 
  parents2  <- as.matrix(training_df[, grep('^Gcov', names(training_df))])/100  # parents in group
  
  parents1   <- parents1[parents1 > 0 ]                             # eliminate zeros 0 so that values can be fit with neg_binom 
  parents2  <- parents2[parents2 > 0 ]                             # eliminate zeros so that values can be fit with neg_binom 
  
  spp_list  <- factor( str_extract( colnames(parents1), '[A-Z]+$'))        # get species names 

  Nspp      <- length(spp_list)                                   # number of parent species
  
  species_name   <- unique(training_df$species )
  
  spp       <- as.numeric( spp_list [spp_list == species_name ] ) # assign species number to datalist 
   
  gid       <- as.numeric(training_df$Group)                      # integer id for each plot area
  G         <- length(unique(training_df$Group))                  # number of groups representing exclosure areas
  
  trackid   <- training_df$trackID
  year      <- training_df$year
  quad      <- training_df$quad
  
  #---------hold out/prediction data ------------------------------------------------------------------
  npreds    <- nrow(holding_df)                                   # total predicted observations, modern data
  y_holdout <- holding_df$Y                                       # new recruit at time t + 1, modern data
  
  Chold     <- holding_df[ , clim_vars]                           # climate matrix, modern data
  
  parents1_out  <- as.matrix(holding_df[, grep('^cov',  names(holding_df)) ])/100  # parents in plot 
  parents2_out  <- as.matrix(holding_df[, grep('^Gcov', names(holding_df))])/100  # parents in group
  
  parents1_out   <- parents1_out[parents1_out > 0 ]                             # eliminate zeros 0 so that values can be fit with neg_binom 
  parents2_out  <- parents2_out[parents2_out > 0 ]                             # eliminate zeros so that values can be fit with neg_binom 
  
  gid_out   <- as.numeric(holding_df$Group)                       # group id, modern data
  yid_out   <- as.numeric(as.factor(holding_df$year))             # year id, modern data
  nyrs_out  <- length(unique(holding_df$year))                    # num years, modern data
  treat_out <- as.numeric(factor(holding_df$Treatment))           # treatments
  
  trackid_out   <- holding_df$trackID
  year_out      <- holding_df$year
  quad_out      <- holding_df$quad
  
  return(
    list(
      N = N, Y = Y, Nspp = Nspp, spp = spp, gid = gid, G = G, Yrs = nyrs, yid = yid , Covs = Covs, C = C, parents1 = parents1, parents2 = parents2,
      gid_out = gid_out, npreds = npreds, y_holdout = y_holdout, Chold = Chold, parents1_out = parents1_out, parents2_out = parents2_out, 
      yid_out = yid_out, nyrs_out = nyrs_out, treat_out = treat_out,
      trackid = trackid, trackid_out = trackid_out, year = year, year_out = year_out, quad = quad, quad_out = quad_out,        
      tau_beta = 1
    )
  )
}


make_stan_datalist <- function(vr, data_path, clim_vars, clim_file, ... ) { 

  clim_vars <- sort(clim_vars ) 

  # -- read data files ---------------------------------------------------------------------# 

  clim_covs <- readRDS(file.path(data_path, clim_file))
  
  dfiles <- dir( data_path, pattern = paste0(vr, '.RDS'), full.names = TRUE)
  
  spp_names <- as.character( regmatches( dfiles, m = gregexpr( pattern = '([A-Z]{4})', dfiles )))
  
  dlist <- lapply( dfiles, readRDS)
  
  # -- subset ------------------------------------------------------------------------------#
  all_data <- lapply(dlist, function(x){ subset(x, year > 1926 & !Treatment %in% c('No_shrub', 'No_grass'))} )
  
  all_data <- lapply(all_data, merge, y = clim_covs, by = c('Treatment', 'Period', 'year')) 
  
  # -- make training and holding subsets ----------------------------------------------------# 
  
  training <- lapply( all_data, function(x) { which(x$Period == "Historical" & x$Treatment == 'Control') } ) 
  holding  <- lapply( all_data, function(x) { which(x$Period == "Modern" )  } ) 

  # -- prepare for stan ---------------------------------------------------------------------# 
  fxn_list <- c('growth_dataframe2datalist', 'survival_dataframe2datalist', 'recruitment_dataframe2datalist')
  
  all_data_list <- mapply( FUN = match.fun( fxn_list [ grep(vr, fxn_list) ]), df = all_data, train = training, hold = holding, MoreArgs = list( 'clim_vars' = clim_vars), SIMPLIFY = FALSE)
  
  names(all_data_list) <- spp_names
  
  # ---- output ------------------------------------------------------------------------------# 
  
  saveRDS(all_data_list, file.path( 'data/temp_data/', paste0( vr, '_data_lists_for_stan.RDS')))
  
}


# -- select covariates -------------------------------------------------------------------#
clim_vars <- c( 'P.a.l', 'P.a.0', 'P.w.sp.0', 'P.w.sp.1', 'T.sp.0', 'T.sp.1', 'T.su.0', 'T.su.1')

clim_file <- 'all_clim_covs.RDS'
data_path <- 'data/temp_data'

make_stan_datalist('survival', data_path, clim_vars, clim_file )

make_stan_datalist('growth', data_path, clim_vars, clim_file )

make_stan_datalist('recruitment', data_path, clim_vars, clim_file )
