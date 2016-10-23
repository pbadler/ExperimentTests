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
    a <- rep(0,nyrs)
    b1_mu <- 1
    b1 <- rep(b1_mu, nyrs)
    sig_a <- sig_a2 <- 0.9
    sig_b1 <- sig_b12 <- 0.13
    bg <- c(0.3, rep(0, G - 1))
    w <- w2 <- -0.15
    tau <- 4
    tau_size <- -0.5
    b2 <- rep( 0, Covs)
    w_all <-c(w, 0, 0, 0)
    a2 <- rep(0, nyrs2)
    b12 <- rep(0, nyrs2)
  }else if ( spp == 'HECO'){
    a <- rep(0,nyrs)
    b1_mu <- 0.9
    b1 <- rep(b1_mu, nyrs)
    sig_a <- sig_a2 <- 0.4
    sig_b1 <- sig_b12 <- 0.1
    bg <- c(0.35, rep(0, G - 1))
    w <- w2 <-  -0.18
    tau <- 0.6
    tau_size <- -0.1
    b2 <- rep( 0, Covs)
    w_all <-c(0, w, 0, 0)
    a2 <- rep(0, nyrs2)
    b12 <- rep(0, nyrs2)
  }else if ( spp == 'POSE'){
    a <- rep(0,nyrs)
    b1_mu <- 0.7
    b1 <- rep(b1_mu, nyrs)
    sig_a <- sig_a2 <- 0.4
    sig_b1 <- sig_b12 <- 0.09
    bg <- c(0.4, rep(0, G - 1))
    w <- w2 <-  -0.4
    tau <- 0.8
    tau_size <- -0.1
    b2 <- rep( 0, Covs)
    w_all <-c(0, 0, w, 0)
    a2 <- rep(0, nyrs2)
    b12 <- rep(0, nyrs2)
  }else if ( spp == 'PSSP'){
    a <- rep(0,nyrs)
    b1_mu <- 1
    b1 <- rep(b1_mu, nyrs)
    sig_a <- sig_a2 <- 1.2
    sig_b1 <- sig_b12 <- 0.25
    bg <- c(0.3, rep(0, G - 1))
    w <- w2 <-  -0.5
    tau <- 0.9
    tau_size <- -0.25
    b2 <- rep( 0, Covs)
    w_all <-c(0, 0, 0, w)
    a2 <- rep(0, nyrs2)
    b12 <- rep(0, nyrs2)
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

dfs <- lapply( dir( 'data/temp_data/', '*scaled_growth_dataframe.RDS', full.names = T), readRDS)
spp <- unlist( lapply( dfs, function(x) unique(x$species)) ) 

# run functions---------------------------------------------------------------------# 
init_vals <- list( NA )

init_vals <- mapply(get_init_vals_growth_models, spp = factor(spp), df = dfs, SIMPLIFY = FALSE)

names(init_vals) <- spp 

# save output ----------------------------------------------------------------------#

saveRDS( init_vals, file = file.path('data/temp_data', 'growth_init_vals.RDS'))
