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

source('R/get_all_demographic_data.R')
source('R/climate/make_climate_variables.R')
source('R/climate/prepare_climate_covariates.R')

split_df <- function(df, train, hold){ 
  
  # split by training and holding 
  training_df <- df[train, ]
  holding_df  <- df[hold, ]
  
  training_df <-  training_df [ order(training_df$year), ]
  holding_df  <-  holding_df [ order(holding_df$year), ]
  
  return( list(training_df, holding_df) ) 
} 

df2list <- function(df, covars, vr, type) { 
  
  N         <- nrow(df)                                           # number of data points for training data 
  gid       <- as.numeric(df$gid)                                 # integer id for each plot area   
  G         <- length(unique(df$gid))                             # number of groups representing exclosure areas
  gm        <- model.matrix.lm(~ df$gid)                         # group contrast matrix 
  
  yid       <- df$yid                                             # integer id for each year 
  nyrs      <- length(unique(df$yid))                             # number of years 

  W         <- as.matrix(df[, grep('W.[A-Z]+', names(df)) ])      # crowding matrix 
  Wcovs     <- ncol(W)                                             # number of species in crowding matrix 
  
  C         <- as.matrix(df[, covars])                            # all climate covariates 
  Covs      <- ncol(C)                                            # number of climate covariates 

  trackid   <- df$trackID
  year      <- df$year
  quad      <- as.numeric(factor(df$quad))
  treat     <- as.numeric(factor(df$Treatment))

  if(vr == 'recruitment'){
    Y         <- df$Y
    parents1  <- as.matrix(df[ , grep('^cov\\.[A-Z]+', names(df))])/100
    parents2  <- as.matrix(df[ , grep('^Gcov\\.[A-Z]+', names(df))])/100
    
    species_name   <- unique(df$species )
    spp_list       <- factor( str_extract( colnames(parents1), '[A-Z]+$')) # get species names 
    spp            <- grep(species_name, spp_list)                  # assign species number to datalist 
    Nspp           <- nlevels(spp_list)
    
  }else { 
    Y         <- df$Y
    X         <- df$X

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
  
  # -------- make year by treatment labels ------------------------------------------------------------ 
  df$treat_year_label <- paste(df$Treatment, df$year, sep = '_')
  df$yid              <- as.numeric(factor(df$treat_year_label)) # each year by treatment gets unique id
  
  # set factors for whole dataset 
  df$gid        <- factor(df$Group)
  df$G          <- nlevels(df$gid)
  
  # get the climate covariates 
  covars <- grep( '^[PT]\\.|^(VWC)\\.', names(df)) 

  if(vr == 'growth'){ 
    df$X <- df$logarea.t0
    df$Y <- df$logarea.t1
  }else if( vr == 'survival'){
    df$X <- df$logarea.t0
    df$Y <- df$survives
  }
  

  # --------split into training and holding data and scale climate covariates -------------------------
  out <- split_df(df, train, hold)
  
  training_df   <- out[[1]]
  holding_df    <- out[[2]]

  # --------training data -----------------------------------------------------------------------------
  training_list  <- df2list(training_df, covars, vr = vr, type = '')
  
  #---------hold out/prediction data ------------------------------------------------------------------
  holding_list  <- df2list(holding_df, covars, vr = vr, type = 'hold')
  
  # save dataframe with covariates ---------------------------------------------
  
  out_df           <- rbind(out[[1]], out[[2]])
  
  saveRDS(out_df, paste0( 'data/temp_data/', species, '_scaled_', vr, '_dataframe.RDS'))
  
  #---------full dataset for estimating year effects ---------------------------------------------------
  full_list <- df2list(out_df, covars, vr = vr, type = '2')
  
  if ( vr == 'growth') { 
    
    #--------use survival dataframe for growth cover predictions ------------------------------# 

    survival_df <- readRDS(paste0('data/temp_data/', species , '_scaled_survival_dataframe.RDS' ))
    survival_df <- survival_df[ survival_df$yid %in% c(full_list$yid2), ] # only use years that are in the growth dataframe  
    survival_df <- subset(survival_df, Period == 'Modern')
    
    covars <- grep( '^[PT]\\.', names(survival_df))     
    cover_list  <- df2list(survival_df, covars, vr = vr, type = '3')
    
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

  all_data_list <- mapply( FUN = compile_datalists, df = all_data, train = training, hold = holding, MoreArgs = list(vr = vr), SIMPLIFY = FALSE)
  
  names(all_data_list) <- spp_names
  
  # ---- output ------------------------------------------------------------------------------# 
  
  saveRDS(all_data_list, file.path( 'data/temp_data/', paste0( vr, '_data_lists_for_stan.RDS')))
  
}


# -- select covariates -------------------------------------------------------------------#
clim_vars <- c('VWC.sp.l_layer1', 'VWC.sp.l_layer2', 
               'VWC.sp.0_layer1', 'VWC.sp.0_layer2', 
               'VWC.sp.1_layer2', 'VWC.sp.1_layer2', 
               'T.sp.0', 'T.sp.1')                     

clim_file <- 'all_clim_covs.RDS'
data_path <- 'data/temp_data'

make_stan_datalist('survival', data_path, clim_vars, clim_file )

make_stan_datalist('growth', data_path, clim_vars, clim_file )

make_stan_datalist('recruitment', data_path, clim_vars, clim_file )