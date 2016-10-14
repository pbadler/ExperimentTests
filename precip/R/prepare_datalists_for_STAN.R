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

scale_covs <- function(df , train, hold){ 
  
  covars <- grep( '^[PT]\\.', names(df)) # grab all the climate covariates for scaling 
  
  training_df <- df[train, ]
  holding_df  <- df[hold, ]
  
  c_means <- colMeans(training_df[ , covars], na.rm = TRUE)
  c_sds   <- apply(training_df[ , covars ], 2, FUN = sd, na.rm = TRUE)
  
  training_df[ , covars] <- scale( training_df[ , covars], center = TRUE, scale = TRUE)
  holding_df[  , covars] <- scale( holding_df[  , covars], center = c_means, scale = c_sds )
  
  training_df <-  training_df [ order(training_df$year), ]
  holding_df  <-  holding_df [ order(holding_df$year), ]
  
  return( list( training_df, holding_df) ) 
} 


growth_dataframe2datalist <- function(df, train, hold){ 
  
  # Function simply makes list of data for STAN models  
  
  survival_dlist <- readRDS('data/temp_data/survival_data_lists_for_stan.RDS')  
  
  # -------- make year by treatment labels ------------------------------------------------------------ 
  df$treat_year_label <- paste(df$Treatment, df$year, sep = '_')
  df$treat_year <- as.numeric(factor(df$treat_year_label))
  
  # --------split into training and holding data and scale climate covariates -------------------------
  out <- scale_covs(df, train, hold)

  covars <- grep( '^[PT]\\.', names(out[[1]])) # grab all the covariates for scaling 
  
  training_df <- out[[1]]
  holding_df <- out[[2]]
  
  # --------training data -----------------------------------------------------------------------------
  N         <- nrow(training_df)                                  # number of data points for training data 
  nyrs      <- length(unique(training_df$treat_year))                   # number of years 
  yid       <- as.numeric(as.factor(training_df$treat_year))            # integer id for each year 
  Y         <- training_df$logarea.t1                             # plant size at time t + 1 
  X         <- training_df$logarea.t0                             # plant size at time t   
  
  C         <- as.matrix(training_df[, covars])                   # all climate covariates 
  Covs      <- ncol(C)                                            # number of climate covariates 
  
  W         <- as.matrix(training_df[, grep('W', names(training_df)) ]) # crowding matrix 
  Wcovs    <- ncol(W)                                             # number of species in crowding matrix 
  
  gid       <- as.numeric(training_df$Group)                      # integer id for each plot area   
  G         <- length(unique(training_df$Group))                  # number of groups representing exclosure areas
  
  trackid   <- training_df$trackID
  year      <- training_df$year
  quad      <- as.numeric(factor(training_df$quad))
  treat     <- as.numeric(factor(training_df$Treatment))
  
  #---------hold out/prediction data ------------------------------------------------------------------
  npreds    <- nrow(holding_df)                                   # total predicted observations, modern data 
  y_holdout <- holding_df$logarea.t1                              # plant size at time t, modern data 
  Xhold     <- holding_df$logarea.t0                              # plant size at time t-1, modern data 
  Chold     <- holding_df[ , covars]                              # climate matrix, modern data 
  Whold     <- holding_df[ , grep('W', names(holding_df)) ]       # crowding matrix, modern data 
  gid_out   <- as.numeric(holding_df$Group)                       # group id, modern data 
  yid_out   <- as.numeric(as.factor(holding_df$treat_year))             # year id, modern data
  nyrs_out  <- length(unique(holding_df$treat_year))                    # num years, modern data 
  treat_out <- as.numeric(factor(holding_df$Treatment))           # treatments 
  
  trackid_out   <- holding_df$trackID
  year_out      <- holding_df$year
  quad_out      <- as.numeric(factor(holding_df$quad))
  
  # save simple dataframe with scaled covariates ---------------------------------------------
  
  species_name   <- unique(out[[1]]$species )
  spp_list       <- factor( str_extract( colnames(W), '[A-Z]+$')) # get species names 
  spp            <- grep(species_name, spp_list)                  # assign species number to datalist 
  out_df         <- do.call(rbind, out)

  out_df$Y       <- out_df$logarea.t1
  out_df$X       <- out_df$logarea.t0
  out_df$yid     <- as.numeric(factor(out_df$treat_year))
  out_df$gid     <- as.numeric(factor(out_df$Group))
  
  saveRDS(out_df, paste0( 'data/temp_data/', species_name, '_scaled_growth_dataframe.RDS'))
  
  #---------full dataset for estimating year effects ---------------------------------------------------
  N2 <- nrow(out_df)
  Y2 <- out_df$Y
  X2 <- out_df$X
  yid2 <- out_df$yid
  Yrs2 <- length(unique(out_df$yid))
  gid2 <- out_df$gid
  W2 <- as.matrix(out_df[, grep('W', names(out_df)) ])
  
  #--------use survival dataframe for cover predictions ------------------------------------------------------------# 

  survival_df <- survival_dlist[[spp]]                      # get data frame for the right species
  
  N3        <- survival_df$npreds                           # total predictions
  X3        <- survival_df$Xhold                            # plant size at time t-1 
  C3        <- survival_df$Chold                            # climate matrix 
  W3        <- survival_df$Whold                            # crowding matrix 
  gid3      <- survival_df$gid_out                          # group id 
  yid3      <- survival_df$yid_out                          # year id
  nyrs3     <- survival_df$nyrs_out                         # num years 
  
  treat3    <- survival_df$treat_out                        # Information for post processing  
  trackid3  <- survival_df$trackid_out
  year3     <- survival_df$year_out
  quad3     <- survival_df$quad_out
  
  return( 
    list(
      N = N, Y = Y, X = X , gid = gid, G = G, Yrs = nyrs, yid = yid , Covs = Covs, C = C, 
      W = W, Whold = Whold,
      Wcovs = Wcovs,
      gid_out = gid_out, npreds = npreds, y_holdout = y_holdout, Xhold = Xhold, Chold = Chold, yid_out = yid_out, nyrs_out = nyrs_out, treat_out = treat_out,
      trackid = trackid, trackid_out = trackid_out, year = year, year_out = year_out, quad = quad, quad_out = quad_out,        
      tau_beta = 1, 
      spp = spp, 
      N2 = N2, Y2 = Y2, X2 = X2, yid2 = yid2, Yrs2 = Yrs2, gid2= gid2, W2 = W2, # for year effects model 
      N3 = N3, X3 = X3, C3 = C3, gid3 = gid3, yid3 = yid3, W3 = W3, # for cover predictions 
      nyrs3 = nyrs3, treat3 = treat3, trackid3 = trackid3, 
      year3 = year3, quad3 = quad3 
    )
  )
}

survival_dataframe2datalist <- function(df, train, hold, covars){

  # Function simply makes list of data for STAN models

  # --------split into training and holding data and scale climate covariates -------------------------
  
  df$treat_year_label <- paste(df$Treatment, df$year, sep = '_')
  df$treat_year <- as.numeric(factor(df$treat_year_label))
  
  out <- scale_covs(df, train, hold)
  
  covars <- grep( '^[PT]\\.', names(out[[1]])) # grab all the covariates for scaling 
  
  training_df <- out[[1]]
  holding_df <- out[[2]]

  # --------training data -----------------------------------------------------------------------------
  N         <- nrow(training_df)                                  # number of data points for training data
  nyrs      <- length(unique(training_df$treat_year))                   # number of years
  yid       <- training_df$treat_year                             # integer id for each year
  Y         <- training_df$survives                               # plant survival at time t + 1   
  X         <- training_df$logarea                                # plant size at time t

  C         <- as.matrix(training_df[, covars])                   # all climate covariates
  Covs      <- ncol(C)                                            # number of climate covariates

  W         <- as.matrix(training_df[, grep('W', names(training_df)) ]) # crowding matrix
  Wcovs     <- ncol(W)                                            # number of crowding effects
  Nspp      <- Wcovs/2
  
  gid       <- as.numeric(training_df$Group)                      # integer id for each plot area
  G         <- length(unique(training_df$Group))                  # number of groups representing exclosure areas

  trackid   <- training_df$trackID
  year      <- training_df$year
  quad      <- as.numeric(factor(training_df$quad))
  treat     <- as.numeric(factor(training_df$Treatment))
  
  #---------hold out/prediction data ------------------------------------------------------------------
  species_name   <- unique(out[[1]]$species )
  spp_list  <- factor( str_extract( colnames(W), '[A-Z]+$'))      # get species names 
  spp       <- grep(species_name, spp_list)                       # assign species number to datalist 
  
  npreds    <- nrow(holding_df)                                   # total predicted observations, modern data
  y_holdout <- holding_df$survives                                # plant survival at time t + 1, modern data
  Xhold     <- holding_df$logarea                                 # plant size at time t, modern data
  Chold     <- holding_df[ , covars]                              # climate matrix, modern data
  Whold     <- holding_df[ , grep('W', names(holding_df)) ]       # crowding matrix, modern data
  gid_out   <- as.numeric(holding_df$Group)                       # group id, modern data
  yid_out   <- as.numeric(as.factor(holding_df$treat_year))             # year id, modern data
  nyrs_out  <- length(unique(holding_df$treat_year))                    # num years, modern data
  treat_out <- as.numeric(factor(holding_df$Treatment))           # treatments

  trackid_out   <- holding_df$trackID
  year_out      <- holding_df$year
  quad_out      <- as.numeric(factor(holding_df$quad))
  
  out_df <-  do.call(rbind, out)
  out_df$Y <- out_df$survives
  out_df$X <- out_df$logarea
  out_df$yid <- as.numeric(factor(out_df$treat_year))
  out_df$gid <- as.numeric(factor(out_df$Group))
  out_df$treat_year_label <- out_df$treat_year_label
  out_df$year <- out_df$year

  saveRDS(out_df, paste0( 'data/temp_data/', species_name, '_scaled_survival_dataframe.RDS'))

  #---------full dataset for estimating year effects ---------------------------------------------------
  
  N2 <- nrow(out_df)
  Y2 <- out_df$Y  
  X2 <- out_df$X
  yid2 <- out_df$yid
  Yrs2 <- length(unique(out_df$yid))
  gid2 <- out_df$gid
  W2 <- as.matrix(out_df[, grep('W', names(out_df)) ])
  C2 <- as.matrix(out_df[, covars])
  treat2 <- as.numeric( factor( out_df$Treatment))
  trackid2 <- out_df$trackID
  year2 <- out_df$year
  quad2 <- as.numeric( factor(out_df$quad) )
  
  return(
    list(
      N = N, Y = Y, X = X , gid = gid, G = G, Yrs = nyrs, yid = yid, 
      Wcovs = Wcovs, W = W,  # competition covariates 
      Covs = Covs, C = C, # climate covariates 
      tau_beta = 1,
      npreds = npreds, y_holdout = y_holdout,  Xhold = Xhold, gid_out = gid_out, yid_out = yid_out, nyrs_out = nyrs_out,
      Whold = Whold, 
      Chold = Chold,
      treat_out = treat_out,
      trackid = trackid, 
      trackid_out = trackid_out,               
      year = year, 
      year_out = year_out,
      quad = quad, 
      quad_out = quad_out,
      spp = spp, 
      Nspp = Nspp, 
      N2 = N2, Y2 = Y2, X2 = X2, yid2 = yid2, Yrs2 = Yrs2, gid2= gid2, W2 = W2, C2 = C2, # for year effects model
      treat2 = treat2, year2 = year2, trackid2 = trackid2, quad2 = quad2
    )
  )
}



recruitment_dataframe2datalist <- function(df, train, hold){
  
  # Function simply makes list of data for STAN models
  
  df$treat_year_label <- paste(df$Treatment, df$year, sep = '_')
  df$treat_year <- as.numeric(factor(df$treat_year_label))
  
  # --------split into training and holding data and scale climate covariates -------------------------
  clim_vars <- names(df)[ grep('^[PT]\\.', names(df)) ] 
  
  out <- scale_covs(df, train, hold)
  
  training_df <- out[[1]]
  holding_df <- out[[2]]

  # --------training data -----------------------------------------------------------------------------
  N         <- nrow(training_df)                                  # number of data points for training data
  
  nyrs      <- length(unique(training_df$treat_year))                   # number of years
  yid       <- as.numeric(as.factor(training_df$treat_year))            # integer id for each year
  Y         <- training_df$Y                                      # new recruit at time t + 1   
  
  C         <- as.matrix(training_df[, clim_vars])                # all climate covariates
  Covs      <- ncol(C)                                            # number of climate covariates
  
  parents1  <- as.matrix(training_df[, grep('^cov', names(training_df)) ])/100  # parents in plot 
  parents2  <- as.matrix(training_df[, grep('^Gcov', names(training_df))])/100  # parents in group
  
  parents2[which(parents2 == 0) ] <- min(parents2[ which(parents2 > 0) ])
  
  spp_list  <- factor( str_extract( colnames(parents1), '[A-Z]+$'))        # get species names 
  Nspp      <- length(spp_list)                                   # number of parent species
  species_name   <- unique(training_df$species )
  spp       <- as.numeric( spp_list [spp_list == species_name ] ) # assign species number to datalist 
   
  gid       <- as.numeric(training_df$Group)                      # integer id for each plot area
  G         <- length(unique(training_df$Group))                  # number of groups representing exclosure areas
  
  trackid   <- training_df$trackID
  year      <- training_df$year
  quad      <- as.numeric( factor( training_df$quad ) ) 
  
  #---------hold out/prediction data ------------------------------------------------------------------
  npreds    <- nrow(holding_df)                                   # total predicted observations, modern data
  y_holdout <- holding_df$Y                                       # new recruit at time t + 1, modern data
  
  Chold     <- holding_df[ , clim_vars]                           # climate matrix, modern data
  
  parents1_out  <- as.matrix(holding_df[, grep('^cov',  names(holding_df)) ])/100  # parents in plot 
  parents2_out  <- as.matrix(holding_df[, grep('^Gcov', names(holding_df))])/100  # parents in group

  parents2_out[which(parents2_out == 0) ] <- min(parents2_out[ which(parents2_out > 0) ])
    
  gid_out   <- as.numeric(holding_df$Group)                       # group id, modern data
  yid_out   <- as.numeric(as.factor(holding_df$treat_year))             # year id, modern data
  nyrs_out  <- length(unique(holding_df$treat_year))                    # num years, modern data
  treat_out <- as.numeric(factor(holding_df$Treatment))           # treatments
  
  trackid_out   <- holding_df$trackID
  year_out      <- holding_df$year
  quad_out      <- as.numeric( factor(holding_df$quad) ) 
  
  out_df <-  do.call(rbind, out)
  out_df$Y <- out_df$Y
  out_df$yid <- as.numeric(factor(out_df$treat_year))
  out_df$gid <- as.numeric(factor(out_df$Group))
  
  saveRDS(out_df, paste0( 'data/temp_data/', species_name, '_scaled_recruitment_dataframe.RDS'))
  
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
  
  # -- subset -------------------------------------------------------------------------------#
  all_data <- lapply(dlist, function(x){ subset(x, year > 1926 & !Treatment %in% c('No_shrub', 'No_grass'))} )
  
  all_data <- lapply(all_data, merge, y = clim_covs[, c('Treatment', 'Period', 'year', clim_vars)], by = c('Treatment', 'Period', 'year')) 
  
  # -- make size by interaction effects -------------------------------------------------------------# 
  
  if(vr == 'survival'){ 
    all_data <- lapply( all_data, function(x) { names(x)[names(x) == 'logarea'] <- 'logarea.t0'; x } )
  } 
  
  if( vr != 'recruitment') { 
    all_data <- 
      lapply( all_data, function( x ) {  
      ifx <- x[, clim_vars]*x[, 'logarea.t0']   ##### climate by size interactions 
      names(ifx ) <- paste0(clim_vars , ':', 'logarea.t0')
      cbind(x, ifx)
    })
    
    all_data <- 
      lapply( all_data, function( x ) {  
        ifx <- x[, grep('W.', names(x))]*x[, 'logarea.t0']   ##### competition by size interactions 
        names(ifx ) <- paste0(names(ifx)[grep('W', names(ifx))] , ':', 'logarea.t0')
        cbind(x, ifx)
    })
  }  
  
  # -- make training and holding subsets ----------------------------------------------------# 
  
  training <- lapply( all_data, function(x) { which(x$Period == "Historical" & x$Treatment == 'Control') } ) 
  holding  <- lapply( all_data, function(x) { which(x$Period == "Modern" )  } ) 

  # -- prepare for stan ---------------------------------------------------------------------# 
  fxn_list <- c('growth_dataframe2datalist', 'survival_dataframe2datalist', 'recruitment_dataframe2datalist')
  
  all_data_list <- mapply( FUN = match.fun( fxn_list [ grep(vr, fxn_list) ]), df = all_data, train = training, hold = holding, SIMPLIFY = FALSE)
  
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
