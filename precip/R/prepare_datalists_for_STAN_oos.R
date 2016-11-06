######################################################################################
#
# Make STAN datalist  
#
#####################################################################################

#rm(list = ls() )

library(stringr)

split_df <- function(df, train, hold){ 
  
  # split by training and holding 
  training_df <- df[train, ]
  holding_df  <- df[hold, ]
  
  training_df <-  training_df [ order(training_df$year), ]
  holding_df  <-  holding_df [ order(holding_df$year), ]
  
  return( list(training_df, holding_df) ) 
} 

scale_covariates <- function( datalist, vr ) { 
  
  datalist$C      <- scale(datalist$C)
  Ccenter         <- attr(datalist$C, 'scaled:center')
  Cscale          <- attr(datalist$C, 'scaled:scale')

  datalist$Chold  <- scale(datalist$Chold, Ccenter, Cscale )

  datalist$W      <- scale(datalist$W)
  Wcenter         <- attr(datalist$W, 'scaled:center')
  Wscale          <- attr(datalist$W, 'scaled:scale')
  
  datalist$Whold  <- scale(datalist$Whold, Wcenter, Wscale )

  if(!is.null( datalist$X)){
      
    X               <- scale(datalist$X)
    datalist$X      <- as.numeric(X)
    Xcenter         <- attr(X, 'scaled:center')
    Xscale          <- attr(X, 'scaled:scale')
  
    datalist$Xhold  <- as.numeric(scale(datalist$Xhold, Xcenter, Xscale))
  }
    
  if( vr == 'growth'){ 
    Y  <- scale(datalist$Y)
    datalist$Y  <- as.numeric(Y)
    Ycenter     <- attr(Y, 'scaled:center')
    Yscale      <- attr(Y, 'scaled:scale') 
    
    datalist$Yhold <- as.numeric(scale(datalist$Yhold, Ycenter, Yscale))

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
  
  yid       <- as.numeric(as.factor( df$yid))                                            # integer id for each year 
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
    
    species_name   <- as.character( unique(df$species ))
    spp_list       <- factor( str_extract( colnames(parents1), '[A-Z]+$')) # get species names 
    spp            <- grep(species_name, spp_list)                  # assign species number to datalist 
    Nspp           <- nlevels(spp_list)
    
  }else { 
    Y         <- df$Y
    X         <- df$X

    species_name   <- as.character( unique(df$species) ) 
    spp_list       <- factor( str_extract( colnames(W), '[A-Z]+$')) # get species names 
    spp            <- grep( species_name, spp_list)                  # assign species number to datalist 
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
  
  out_list <-   c(training_list, holding_list)
  
  out_list <-  scale_covariates(out_list, vr = vr )

  return(out_list)

} 
  
  
make_stan_datalist <- function(vr, data_path, clim_vars, clim_file, year_oos ) { 

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
  
  training <- lapply( all_data, function(x) { which(x$Period == "Historical" & x$Treatment == 'Control' & !x$year %in% year_oos) } ) 
  holding  <- lapply( all_data, function(x) { which(x$year %in% year_oos )  } ) 

  # -- prepare for stan ---------------------------------------------------------------------# 

  all_data_list <- mapply( FUN = compile_datalists, df = all_data, train = training, hold = holding, MoreArgs = list(vr = vr), SIMPLIFY = FALSE)
  
  names(all_data_list) <- spp_names
  
  # ---- output ------------------------------------------------------------------------------# 
  
  return(all_data_list)
  
}


# -- select covariates -------------------------------------------------------------------#
clim_vars <- c('VWC.sp.l', 
               'VWC.sp.0', 
               'VWC.sp.1',
               'VWC.su.0', 
               'VWC.su.l', 
               'T.sp.0', 
               'T.sp.1', 
               'T.sp.l')                     

clim_file <- 'all_clim_covs.RDS'
#data_path <- 'data/temp_data'

# make_stan_datalist('survival', data_path, clim_vars, clim_file )
# 
# make_stan_datalist('growth', data_path, clim_vars, clim_file )
# 
# make_stan_datalist('recruitment', data_path, clim_vars, clim_file )