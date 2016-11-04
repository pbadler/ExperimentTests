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

get_init_vals_survival_models <- function( spp, df, ... ) {
  
  df_old <- subset(df, Period == 'Historical')
  nyrs2 <- length(unique(df$treat_year))
  nyrs <- length(unique(df_old$yid))
  G <- length(unique(df$gid))
  Covs <- length(df[, grep('^[T]\\.|^(VWC)\\.', names(df))])
  
  if(spp == 'ARTR'){
    a_raw <- rep(0,nyrs)
    b1_mu <- 1
    b1 <- rep(b1_mu, nyrs)
    b1_raw <- rep(0, nyrs)
    sig_a <- 0.5
    sig_b1 <- 0.18
    bg <- c(-2.3, rep(0, G - 1))
    w <- -1.5
    b2 <- rep(0, Covs)
    w_all <-c(w, 0, 0, 0)
  }else if ( spp == 'HECO'){
    a_raw <- rep(0,nyrs)
    b1_mu <- 1
    b1 <- rep(b1_mu, nyrs)
    b1_raw <- rep(0, nyrs)
    sig_a <- 0.5
    sig_b1 <- 0.18
    bg <- c(-2.3, rep(0, G - 1))
    w <- -1.5
    b2 <- rep(0, Covs)
    w_all <-c(0, w, 0, 0)
  }else if ( spp == 'POSE'){
    a_raw <- rep(0,nyrs)
    b1_mu <- 1
    b1 <- rep(b1_mu, nyrs)
    b1_raw <- rep(0, nyrs)
    sig_a <- 0.5
    sig_b1 <- 0.18
    bg <- c(-2.3, rep(0, G - 1))
    w <- -1.5
    b2 <- rep(0, Covs)
    w_all <-c(0, 0, w, 0)
  }else if ( spp == 'PSSP'){
    a_raw <- rep(0,nyrs)
    b1_mu <- 1
    b1 <- rep(b1_mu, nyrs)
    b1_raw <- rep(0, nyrs)
    sig_a <- 0.5
    sig_b1 <- 0.18
    bg <- bg2 <- c(-2.3, rep(0, G - 1))
    w <- -1.5
    b2 <- rep(0, Covs)
    w_all <-c(0, 0, 0, w)
  }
  
  rm(df, df_old)
  spp <- as.numeric(factor(spp))
  
  init_vals <-  lapply( ls(), function(x) eval(parse( text = x) ))  # collect inits 
  names( init_vals) <- ls()[-which(ls() == 'init_vals')]
  
  init_vals$w <- init_vals$w_all
  
  return(init_vals)
}

# input files ----------------------------------------------------------------------#

dfs <- lapply( dir( 'data/temp_data/', '*_survival_cleaned_dataframe.RDS', full.names = T), readRDS)
spp <- unlist( lapply( dfs, function(x) unique(x$species)) ) 

# run functions---------------------------------------------------------------------# 
init_vals <- list( NA )

init_vals <- mapply(get_init_vals_survival_models, spp = factor(spp), df = dfs, SIMPLIFY = FALSE)

names(init_vals) <- spp 

# save output ----------------------------------------------------------------------#

saveRDS( init_vals, file = file.path('data/temp_data', 'survival_init_vals.RDS'))

