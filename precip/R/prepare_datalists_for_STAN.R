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

scale_covs <- function(df, train, hold){ 
  
  # -------- make year by treatment labels ------------------------------------------------------------ 
  df$treat_year_label <- paste(df$Treatment, df$year, sep = '_')
  df$yid              <- as.numeric(factor(df$treat_year_label)) # each year by treatment gets unique id
  
  # set factors for whole dataset 
  df$gid        <- factor(df$Group)
  df$G          <- nlevels(df$gid)

  # get the climate covariates 
  covars <- grep( '^[PT]\\.', names(df)) 
  
  training_df <- df[train, ]
  holding_df  <- df[hold, ]
  
  c_means <- colMeans(training_df[ , covars], na.rm = TRUE)
  c_sds   <- apply(training_df[ , covars ], 2, FUN = sd, na.rm = TRUE)
  
  training_df[ , covars] <- scale( training_df[ , covars], center = TRUE, scale = TRUE)
  holding_df[  , covars] <- scale( holding_df[  , covars], center = c_means, scale = c_sds )
  
  training_df <-  training_df [ order(training_df$year), ]
  holding_df  <-  holding_df [ order(holding_df$year), ]
  
  return( list(training_df, holding_df, covars) ) 
} 

df2list <- function(df, covars, vr, type) { 
  
  N         <- nrow(df)                                           # number of data points for training data 
  gid       <- as.numeric(df$gid)                                 # integer id for each plot area   
  G         <- length(unique(df$gid))                             # number of groups representing exclosure areas
  gm        <- model.matrix.lm(~ df$gid)                         # group contrast matrix 
  
  yid       <- df$yid                                             # integer id for each year 
  nyrs      <- length(unique(df$yid))                             # number of years 

  W         <- as.matrix(df[, grep('W.[A-Z]+', names(df)) ])      # crowding matrix 
  Wcovs    <- ncol(W)                                             # number of species in crowding matrix 
  
  C         <- as.matrix(df[, covars])                            # all climate covariates 
  Covs      <- ncol(C)                                            # number of climate covariates 

  trackid   <- df$trackID
  year      <- df$year
  quad      <- as.numeric(factor(df$quad))
  treat     <- as.numeric(factor(df$Treatment))

  if(vr == 'recruitment'){
    Y         <- df$Y
    parents1  <- df$parents1
    parents2  <- df$parents2 
    
    species_name   <- unique(df$species )
    spp_list       <- factor( str_extract( colnames(parents1), '[A-Z]+$')) # get species names 
    spp            <- grep(species_name, spp_list)                  # assign species number to datalist 
    
  }else if( vr == 'survival'){ 
    Y         <- df$survives
    X         <- df$logarea.t0
    W         <- as.matrix(df[, grep('W\\.[A-Z]+', names(df)) ])      # crowding matrix 
    Wcovs     <- ncol(W)                                               # number of species in crowding matrix
    
    species_name   <- unique(df$species)
    spp_list       <- factor( str_extract( colnames(W), '[A-Z]+$')) # get species names 
    spp            <- grep(species_name, spp_list)                  # assign species number to datalist 
    
  }else if( vr == 'growth'){ 
    Y         <- df$logarea.t1
    X         <- df$logarea.t0
    W         <- as.matrix(df[, grep('W.[A-Z]+', names(df)) ])      # crowding matrix 
    Wcovs     <- ncol(W)                                             # number of species in crowding matrix 
    
    species_name   <- unique(df$species)
    spp_list       <- factor( str_extract( colnames(W), '[A-Z]+$')) # get species names 
    spp            <- grep(species_name, spp_list)                  # assign species number to datalist 
    
  }
  
  out <- lapply( ls(), function(x) eval(parse(text = x)))
  lnames <- ls()[- grep('out', ls()) ]
  names(out) <- paste0(lnames, type)  
  out <- out[- grep('df', names(out))]
  
  return(out)
}

compile_datalists <- function( df, train, hold, vr ) { 

  species <- as.character(  unique(df$species) )
  
  # Function simply makes list of data for STAN models  
  
  # --------split into training and holding data and scale climate covariates -------------------------
  out <- scale_covs(df, train, hold)
  
  training_df   <- out[[1]]
  holding_df    <- out[[2]]
  covars        <- out[[3]]
  
  # --------training data -----------------------------------------------------------------------------
  training_list  <- df2list(training_df, covars, vr = vr, type = '')
  
  #---------hold out/prediction data ------------------------------------------------------------------
  holding_list  <- df2list(holding_df, covars, vr = vr, type = 'hold')
  
  # save simple dataframe with scaled covariates ---------------------------------------------
  
  out_df           <- rbind(out[[1]], out[[2]])
  if(vr == 'survival'){ 
    out_df$Y       <- out_df$survives
  }else { 
    out_df$Y       <- out_df$logarea.t1
  }
  out_df$X         <- out_df$logarea.t0
  
  saveRDS(out_df, paste0( 'data/temp_data/', species, '_scaled_', vr, '_dataframe.RDS'))
  
  #---------full dataset for estimating year effects ---------------------------------------------------
  full_list <- df2list(out_df, covars, vr = vr, type = '2')
  
  if ( vr == 'growth') { 
    #--------use survival dataframe for growth cover predictions ------------------------------# 
    
    survival_dlist <- readRDS('data/temp_data/survival_data_lists_for_stan.RDS')  
    
    survival_df <- survival_dlist[[species]]                      # get data frame for the right species
    
    cover_list <- list()
    cover_list$N3        <- survival_df$npreds                           # total predictions
    cover_list$X3        <- survival_df$Xhold                            # plant size at time t-1 
    cover_list$C3        <- survival_df$Chold                            # climate matrix 
    cover_list$W3        <- survival_df$Whold                            # crowding matrix 
    cover_list$gid3      <- survival_df$gid_out                          # group id 
    cover_list$yid3      <- survival_df$yid_out                          # year id
    cover_list$nyrs3     <- survival_df$nyrs_out                         # num years 
    
    cover_list$treat3    <- survival_df$treat_out                        # Information for post processing  
    cover_list$trackid3  <- survival_df$trackid_out
    cover_list$year3     <- survival_df$year_out
    cover_list$quad3     <- survival_df$quad_out
    
    out_list <-   c(training_list, holding_list, full_list, cover_list)
    
  } else { 
    
  out_list <-   c(training_list, holding_list, full_list)
  
  }
  
  return(out_list)

} 
  
  
make_stan_datalist <- function(vr, data_path, clim_vars, clim_file ) { 

  clim_vars <- sort(clim_vars ) 

  # -- read data files ---------------------------------------------------------------------# 

  clim_covs <- readRDS(file.path(data_path, clim_file))
  
  dfiles <- dir( data_path, pattern = paste0(vr, '.RDS'), full.names = TRUE)
  
  spp_names <- as.character( regmatches( dfiles, m = gregexpr( pattern = '([A-Z]{4})', dfiles )))
  dlist <- lapply( dfiles, readRDS)
  
  # -- get indeces for training and holding data ----------------------------------------------------#
  all_data <- lapply(dlist, function(x){ subset(x, !Treatment %in% c('No_shrub', 'No_grass'))} )
  
  all_data <- lapply(all_data, merge, y = clim_covs[, c('Treatment', 'Period', 'year', clim_vars)], by = c('Treatment', 'Period', 'year')) 
  
  if(vr == 'survival'){ 
    all_data <- lapply( all_data, function(x) { names(x)[names(x) == 'logarea'] <- 'logarea.t0'; x } )
  } 
  
  # -- make size by interaction effects -------------------------------------------------------------# 
  
  if( vr != 'recruitment') {
    all_data <-
      lapply( all_data, function( x ) {
          ifx <- x[, clim_vars]*x[, 'logarea.t0']   ##### climate by size interactions
          names(ifx ) <- paste0(clim_vars , ':', 'logarea.t0')
          cbind(x, ifx)
        }
      )
  #   # all_data <-
  #   #   lapply( all_data, function( x ) {
  #   #     ifx <- x[, grep('W.', names(x))]*x[, 'logarea.t0']   ##### competition by size interactions
  #   #     names(ifx ) <- paste0(names(ifx)[grep('W.[A-Z]+', names(ifx))] , ':', 'logarea.t0')
  #   #     cbind(x, ifx)
  #   #   })
  }

  # -- make training and holding indeces ----------------------------------------------------# 
  
  training <- lapply( all_data, function(x) { which(x$Period == "Historical" & x$Treatment == 'Control') } ) 
  holding  <- lapply( all_data, function(x) { which(x$Period == "Modern" )  } ) 

  # -- prepare for stan ---------------------------------------------------------------------# 
  #fxn_list <- c('growth_dataframe2datalist', 'survival_dataframe2datalist', 'recruitment_dataframe2datalist')
  
  all_data_list <- mapply( FUN = compile_datalists, df = all_data, train = training, hold = holding, MoreArgs = list(vr = vr), SIMPLIFY = FALSE)
  
  names(all_data_list) <- spp_names
  
  # ---- output ------------------------------------------------------------------------------# 
  
  saveRDS(all_data_list, file.path( 'data/temp_data/', paste0( vr, '_data_lists_for_stan.RDS')))
  
}


# -- select covariates -------------------------------------------------------------------#
clim_vars <- c( 'P.f.w.sp.l', 'P.f.w.sp.0', 'P.f.w.sp.1', 'P.su.0', 'P.su.1', 'T.sp.0', 'T.sp.1')

clim_file <- 'all_clim_covs.RDS'
data_path <- 'data/temp_data'

make_stan_datalist('survival', data_path, clim_vars, clim_file )

make_stan_datalist('growth', data_path, clim_vars, clim_file )

make_stan_datalist('recruitment', data_path, clim_vars, clim_file )

