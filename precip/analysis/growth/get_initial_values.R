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

get_init_vals_growth_models <- function( spp, df, ... ) {
  
  df_old <- subset(df, Period == 'Historical')
  nyrs2 <- length(unique(df$treat_year))
  nyrs <- length(unique(df_old$yid))
  G <- length(unique(df$gid))
  Covs <- length(df[, grep('^[PT]\\.', names(df))])
  
  if(spp == 'ARTR'){
    sigma <- sigma2 <- 0.75
    a_raw <- rep(0,nyrs)
    b1_mu <- 1
    b1_mu2 <- 0.9
    b1 <- rep(b1_mu, nyrs)
    b1_raw <- rep(0, nyrs)
    sig_a <- sig_a2 <- 1
    sig_b1 <- sig_b12 <- 0.18
    bg <- bg2 <- c(0.7, rep(0, G - 1))
    w <- w2 <- -0.15
    b2 <- rep( 0, Covs)
    w_all <-c(w, 0, 0, 0)
    a_raw2 <- rep(0,nyrs2)
    b1_raw2 <- rep(0, nyrs2)
  }else if ( spp == 'HECO'){
    sigma <- sigma2 <- 0.8
    a_raw <- rep(0,nyrs)
    b1_mu <- 0.9
    b1_mu2 <- 0.8
    b1 <- rep(b1_mu, nyrs)
    b1_raw <- rep(0, nyrs)
    sig_a <- sig_a2 <- 0.4
    sig_b1 <- sig_b12 <- 0.1
    bg <- bg2 <- c(0.35, rep(0, G - 1))
    w <- w2 <-  -0.05
    b2 <- rep( 0, Covs)
    w_all <-c(0, w, 0, 0)
    a_raw2 <- rep(0,nyrs2)
    b1_raw2 <- rep(0, nyrs2)
  }else if ( spp == 'POSE'){
    sigma <- sigma2 <- 1.02
    a_raw <- rep(0,nyrs)
    b1_mu <- 0.35
    b1_mu2 <- 0.67
    b1 <- rep(b1_mu, nyrs)
    b1_raw <- rep(0, nyrs)
    sig_a <- sig_a2 <- 0.4
    sig_b1 <- sig_b12 <- 0.09
    bg <- bg2 <- c(0.7, rep(0, G - 1))
    w <- w2 <-  -0.4
    b2 <- rep( 0, Covs)
    w_all <-c(0, 0, w, 0)
    a_raw2 <- rep(0,nyrs2)
    b1_raw2 <- rep(0, nyrs2)
  }else if ( spp == 'PSSP'){
    sigma <- sigma2 <- 0.85
    a_raw <- rep(0,nyrs)
    b1_mu <- 1
    b1_mu2 <- 0.8
    b1 <- rep(b1_mu, nyrs)
    b1_raw <- rep(0, nyrs)
    sig_a <- sig_a2 <- 0.35
    sig_b1 <- sig_b12 <- 0.1
    bg <- bg2 <- c(0.3, rep(0, G - 1))
    w <- w2 <-  -0.5
    b2 <- rep( 0, Covs)
    w_all <-c(0, 0, 0, w)
    a_raw2 <- rep(0,nyrs2)
    b1_raw2 <- rep(0, nyrs2)
  }
  
  rm(df, df_old)
  spp <- as.numeric(spp)
  
  init_vals <-  lapply( ls(), function(x) eval(parse( text = x) ))  # collect inits 
  names( init_vals) <- ls()[-which(ls() == 'init_vals')]

  inits <- rep( list(init_vals), 3 ) 
  
  inits[[3]]$w <- init_vals$w_all
  inits[[3]]$w2 <- init_vals$w_all
  
  return(inits)
}

# input files ----------------------------------------------------------------------#

dfs <- lapply( dir( 'data/temp_data/', '*scaled_growth_dataframe.RDS', full.names = T), readRDS)
spp <- unlist( lapply( dfs, function(x) unique(x$species)) ) 

# run functions---------------------------------------------------------------------# 
init_vals <- list( NA )

init_vals <- mapply(get_init_vals_growth_models, spp = factor(spp), df = dfs, SIMPLIFY = FALSE)

names(init_vals) <- spp 

# save output ----------------------------------------------------------------------#

saveRDS( init_vals, file = file.path('data/temp_data', 'growth_init_vals.RDS'))
