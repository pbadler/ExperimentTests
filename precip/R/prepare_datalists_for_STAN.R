######################################################################################
#
# Make STAN datalist  
#
#####################################################################################

rm(list = ls() )

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

detachAllPackages()

# source('R/ExtractData_3Runs.R')
# source('R/aggregate_spot_VWC.R')
# source('R/soilMoistureTreatmentEffects.R')
# source('R/climate/aggregate_VWC_data.R')
# source('R/get_all_demographic_data.R')
# source('R/climate/make_climate_variables.R')
# source('R/climate/prepare_climate_covariates.R')

library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(stringr)

split_df <- function(df, train, hold){ 
  
  # split by training and holding 
  training_df <- df[train, ]
  holding_df  <- df[hold, ]
  
  training_df <-  training_df [ order(training_df$year), ]
  holding_df  <-  holding_df [ order(holding_df$year), ]
  
  return( list(training_df, holding_df) ) 
} 

scale_covariates <- function( datalist) { 
  
  datalist$C      <- scale(datalist$C)
  Ccenter         <- attr(datalist$C, 'scaled:center')
  Cscale          <- attr(datalist$C, 'scaled:scale')

  datalist$Chold  <- scale(datalist$Chold, Ccenter, Cscale )
  datalist$C2     <- scale(datalist$C2, Ccenter, Cscale ) 

  datalist$W      <- scale(datalist$W)
  Wcenter         <- attr(datalist$W, 'scaled:center')
  Wscale          <- attr(datalist$W, 'scaled:scale')
  
  datalist$Whold  <- scale(datalist$Whold, Wcenter, Wscale )
  datalist$W2     <- scale(datalist$W2, Wcenter, Wscale ) 
  
  datalist$Ccenter <- Ccenter 
  datalist$Cscale  <- Cscale 
  datalist$Wcenter <- Wcenter
  datalist$Wscale  <- Wscale 
  
  if(!is.null( datalist$X)){
    
    # ifx <- datalist$C*datalist$X
    # 
    # names(ifx) <- paste0(colnames(datalist$C), 'x', 'logarea.t0')
    # datalist$C <- cbind(datalist$C, ifx)
    #       
    # ifx2 <- datalist$C2*datalist$X2
    # names(ifx2) <- paste0(colnames(datalist$C2), 'x', 'logarea.t0')
    # datalist$C2 <- cbind(datalist$C2, ifx2)
    # 
    # ifxh <- datalist$Chold*datalist$Xhold
    # names(ifxh) <- paste0(colnames(datalist$Chold), 'x', 'logarea.t0')
    # datalist$Chold <- cbind(datalist$Chold, ifxh)
    
    X               <- scale(datalist$X)
    datalist$X      <- as.numeric(X)
    Xcenter         <- attr(X, 'scaled:center')
    Xscale          <- attr(X, 'scaled:scale')

    datalist$Xhold  <- as.numeric(scale(datalist$Xhold, Xcenter, Xscale))
    datalist$X2     <- as.numeric(scale(datalist$X2, Xcenter, Xscale))

    datalist$Xcenter <- Xcenter
    datalist$Xscale  <- Xscale
    
  }
    
  if(!is.null( datalist$C3)){
    datalist$C3 <- scale(datalist$C3, Ccenter, Cscale)
    datalist$W3 <- scale(datalist$W3, Wcenter, Wscale)
    datalist$X3 <- as.numeric(scale(datalist$X3, Xcenter, Xscale))
    
    # ifx3 <- datalist$C3*datalist$X3
    # names(ifx3) <- paste0(colnames(datalist$C3), 'x', 'logarea.t0')
    # datalist$C3 <- cbind(datalist$C3, ifx3)
    # 
    # Y  <- scale(datalist$Y)
    # datalist$Y  <- as.numeric(Y)
    # Ycenter     <- attr(Y, 'scaled:center')
    # Yscale      <- attr(Y, 'scaled:scale') 
    # 
    # datalist$Yhold <- as.numeric(scale(datalist$Yhold, Ycenter, Yscale))
    # datalist$Y2 <- as.numeric(scale(datalist$Y2, Ycenter, Yscale))
    # 
    # datalist$Ycenter <- Ycenter
    # datalist$Yscale  <- Yscale 
  
  }

  class_list    <-  unlist( lapply( datalist, function(x) class(x)))
   
  datalist <- datalist[ which(class_list != 'character')]  
    
  return( datalist )
}


df2list <- function(df, covars, vr, type) { 
  
  N         <- nrow(df)                                           # number of data points for training data 
  gid       <- as.numeric(df$gid)                                 # integer id for each plot area   
  G         <- length(unique(df$gid))                             # number of groups representing exclosure areas
  gm        <- model.matrix.lm(~ df$gid)                         # group contrast matrix 
  
  yid       <- df$yid                                             # integer id for each year 
  nyrs      <- length(unique(df$yid))                             # number of years 

  W         <- as.matrix(df[, grep('W.[A-Z]+', names(df)) ])      # crowding matrix 
  Wcovs     <- ncol(W)                                            # number of species in crowding matrix 
  
  C         <- as.matrix(df[, covars])                            # all climate covariates 
  Covs      <- ncol(C)                                            # number of climate covariates 

  trackid   <- df$trackID
  year      <- df$year
  quad      <- as.numeric(factor(df$quad))
  treat     <- as.numeric(factor(df$Treatment))

  if( type %in% c('2', 'hold')) { 
    tm        <- model.matrix.lm(~ df$Treatment)[, -1]                  # Treatment matrix, no intercept
    nT        <- ncol(tm)                                         # number of treatments 
  }
  
  if(vr == 'recruitment'){
    Y              <- df$Y
    parents1       <- as.matrix(df[ , grep('^cov\\.[A-Z]+', names(df))])/100
    parents2       <- as.matrix(df[ , grep('^Gcov\\.[A-Z]+', names(df))])/100
    
    species_name   <- unique(df$species )
    spp_list       <- factor( str_extract( colnames(parents1), '[A-Z]+$')) # get species names 
    spp            <- grep(species_name, spp_list)                  # assign species number to datalist 
    Nspp           <- nlevels(spp_list)
    
  }else { 
    Y         <- df$Y
    X         <- df$X
    
    if( type %in% c('2', 'hold')) { 
      tm        <- cbind( tm , tm*X ) # treatment by size interaction 
      nT        <- ncol(tm)
    }
    
    species_name   <- unique(df$species)
    spp_list       <- factor( str_extract( colnames(W), '[A-Z]+$')) # get species names 
    spp            <- grep(species_name, spp_list)                  # assign species number to datalist 
  }    
  
  rm(vr, species_name, spp_list)

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
  
  saveRDS(out_df, paste0( 'data/temp_data/', species,'_', vr, '_cleaned_dataframe.RDS'))
  
  #---------full dataset for estimating year effects ---------------------------------------------------
  full_list <- df2list(out_df, covars, vr = vr, type = '2')
  
  if ( vr == 'growth') { 
    
    #--------use survival dataframe for growth cover predictions ------------------------------# 

    survival_df <- readRDS(paste0('data/temp_data/', species , '_survival_cleaned_dataframe.RDS' ))

    survival_df <- survival_df[ survival_df$yid %in% c(full_list$yid2), ] # only use years that are in the growth dataframe  
    
    #survival_df <- subset(survival_df, Period == 'Modern')
    covars <- grep( '^[PT]\\.|^(VWC)\\.', names(survival_df)) 

    cover_list  <- df2list(survival_df, covars, vr = vr, type = '3')
    
    out_list <-   c(training_list, holding_list, full_list, cover_list)
    
  } else { 
    
  out_list <-   c(training_list, holding_list, full_list)
  
  }
  
  out_list <-  scale_covariates(out_list)

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
  
  all_data <- lapply( all_data, function( x ) x[complete.cases(x), ] ) # remove rows with NAs
  
  if(vr == 'survival'){ 
    all_data <- lapply( all_data, function(x) { names(x)[names(x) == 'logarea'] <- 'logarea.t0'; x } )
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
clim_vars <- c('VWC.sp.l', 
               'VWC.sp.0', 
               'VWC.sp.1',
               'VWC.su.l', 
               'VWC.su.0', 
               'VWC.su.1',
               'VWC.f.l', 
               'VWC.f.0', 
               'VWC.f.1',
               'T.sp.1', 
               'T.sp.0',
               'T.sp.l',
               'T.su.1', 
               'T.su.0', 
               'T.su.l',
               'T.f.1', 
               'T.f.0', 
               'T.f.l')                     

clim_file <- 'all_clim_covs.RDS'
data_path <- 'data/temp_data'

make_stan_datalist('survival', data_path, clim_vars, clim_file )

make_stan_datalist('growth', data_path, clim_vars, clim_file )

make_stan_datalist('recruitment', data_path, clim_vars, clim_file )
