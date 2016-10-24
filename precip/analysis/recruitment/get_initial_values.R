##################################################################################################################
#
# Get init vals for stan models
#
##################################################################################################################

#---model descriptions --------------------------------------------------------------------------------------------
#  
#   m1: null model:                       w/ size  + intra-specific competition                      
#   m2: single species model              w/ size  + intra-specific competition + climate 
#   m3: multi species model               w/ size  + all competition + climate  
# ----------------------------------------------------------------------------------------------------------------- 

rm(list = ls())

library(lme4)
library(parallel)

get_init_vals_recruitment_models <- function( spp, df, ... ) {
  
  df_old <- subset(df, Period == 'Historical')
  nyrs2 <- length(unique(df$treat_year))
  nyrs <- length(unique(df_old$yid))
  G <- length(unique(df$gid))
  Covs <- length(df[, grep('^[PT]\\.', names(df))])
  
  if(spp == 'ARTR'){
    bg  <- c(-2, rep(0, G - 1))
    a_raw <- rep(0,nyrs)
    sig_a <- sig_a2 <- 1.3
    w <- w2 <-  -2
    b2 <- rep( 0, Covs)
    w_all <-c(0, 0, 0, w)
    a_raw2 <- rep(0,nyrs2)
    u <- 0.9
    u2 <- 0.9
    theta <- 1.2
    theta2 <- 1.2
    
    
  }else if ( spp == 'HECO'){
    bg  <- c(-2, rep(0, G - 1))
    a_raw <- rep(0,nyrs)
    sig_a <- sig_a2 <- 1.3
    w <- w2 <-  -2
    b2 <- rep( 0, Covs)
    w_all <-c(0, 0, 0, w)
    a_raw2 <- rep(0,nyrs2)
    u <- 0.9
    u2 <- 0.9
    theta <- 1.2
    theta2 <- 1.2
    
    
  }else if ( spp == 'POSE'){
    bg  <- c(-2, rep(0, G - 1))
    a_raw <- rep(0,nyrs)
    sig_a <- sig_a2 <- 1.3
    w <- w2 <-  -2
    b2 <- rep( 0, Covs)
    w_all <-c(0, 0, 0, w)
    a_raw2 <- rep(0,nyrs2)
    u <- 0.9
    u2 <- 0.9
    theta <- 1.2
    theta2 <- 1.2
    
    
  }else if ( spp == 'PSSP'){
    bg  <- c(-2, rep(0, G - 1))
    a_raw <- rep(0,nyrs)
    sig_a <- sig_a2 <- 1.3
    w <- w2 <-  -2
    b2 <- rep( 0, Covs)
    w_all <-c(0, 0, 0, w)
    a_raw2 <- rep(0,nyrs2)
    u <- 0.9
    u2 <- 0.9
    theta <- 1.2
    theta2 <- 1.2
    
  }
  
  rm(df, df_old)
  
  init_vals <-  lapply( ls(), function(x) eval(parse( text = x) ))  # collect inits 
  names( init_vals) <- ls()[-which(ls() == 'init_vals')]
  
  inits <- rep( list(init_vals), 3 ) 
  
  inits[[3]]$w <- init_vals$w_all_2
  inits[[3]]$w2 <- init_vals$w_all_2
  
  return(inits)
}

# input files ----------------------------------------------------------------------#

dfs <- lapply( dir( 'data/temp_data/', '*scaled_recruitment_dataframe.RDS', full.names = T), readRDS)
spp <- unlist( lapply( dfs, function(x) unique(x$species)) ) 

# run functions---------------------------------------------------------------------# 
init_vals <- list( NA )

init_vals <- mapply(get_init_vals_recruitment_models, spp = factor(spp), df = dfs, SIMPLIFY = FALSE)

names(init_vals) <- spp 

# save output ----------------------------------------------------------------------#

saveRDS( init_vals, file = file.path('data/temp_data', 'recruitment_init_vals.RDS'))


