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
#   m4: year effects model                w/ size  + all competition + year effects 
# ----------------------------------------------------------------------------------------------------------------- 

rm(list = ls())

library(lme4)
library(parallel)

set_init_vals_list <-  function( model, C_names, W_names , nyrs2) {  
  
  init_vals <- as.list( fixef(model)[1] )
  
  init_vals <- c(init_vals, as.numeric(data.frame( VarCorr(model) )$sdcor[ c(1,2)]))
  
  names( init_vals )[1:3] <- c('a_mu', 'sig_a', 'sig_G')
  
  if(init_vals$sig_G < 0.0001) { init_vals$sig_G <- 0.1 } # prevent sig G from being zero in some cases 
  
  b2 <- fixef(model)[C_names]
  w <- fixef(model)[W_names]

  if( length(b2[!is.na(b2)]) > 0 )  { 
    init_vals <- c(init_vals, list( b2 = as.numeric(b2[!is.na(b2)])))
  } 
  if( length(w[!is.na(W_names)]) > 0 )  { 
    init_vals <- c(init_vals, list(w = as.numeric(w[!is.na(w)])))  
  }
  
  # random effects initial values need to be specified for STAN 2.6.0
  nyrs <- length(unique(model@frame$yid) )
  G <- length(unique(model@frame$gid)) 
  
  init_vals$a <- rep(0, nyrs)
  init_vals$gint <- rep(0, G) 
  
  init_vals$theta <- median(model@theta)
  init_vals$u <- 0.9
  
  return( init_vals )

}


get_init_vals_recruitment_models <- function( spp, df, ... ) {
  
  df_old <- subset(df, Period == 'Historical')

  nyrs2 <- length(unique(df$treat_year))
  
  nyrs <- length(unique(df_old$yid))
  G <- length(unique(df$gid))
  Ccovs <- length(df[, grep('^[PT]\\.', names(df))])
  
  if(spp == 'ARTR'){
    a_mu <- a_mu2 <- 1.2
    sig_a <- sig_a2 <- 1.3
    sig_G <- sig_G2 <- 0.5
    w <- w2 <- -0.4
    a <- rep(a_mu, nyrs)
    a2 <- rep(a_mu, nyrs2)
    gint <- gint2 <- rep( 0, G)
    theta <- theta2 <- 0.6
    u <- u2 <- 0.75
    b2 <- c(-0.4, -0.8, 0.7, 1, 0.5, 0, -0.1, -0.9)
    w_all <- w_all_2 <- c(w, -0.1, -0.1, -0.1)
    
  }else if(spp == 'HECO'){
    a_mu <- a_mu2 <- 2.6
    sig_a <- sig_a2 <- 1.2
    sig_G <- sig_G2 <- 0.1
    w <- w2 <-  -1.6
    a <- rep(a_mu, nyrs)
    a2 <- rep(a_mu, nyrs2)
    gint <- gint2 <- rep( 0, G)
    theta <- theta2 <- 1.2
    u <- u2 <-  0.90
    b2 <- c(-0.3, 0.04, 0.3, 0.4, 0.2, 0.5, 0.3, -0.1)
    
    w_all <- w_all_2 <- c(-0.1, w, -0.1, -0.1)
    
 }else if ( spp == 'POSE'){
    a_mu <- a_mu2 <- 3.4
    sig_a <- sig_a2 <- 0.8
    sig_G <- sig_G2 <- 0.3
    w <- w2 <-  -2.1
    a <- rep(a_mu, nyrs)
    a2 <- rep(a_mu, nyrs2)
    gint <- gint2 <- rep( 0, G)
    theta <- theta2 <- 1.23
    u <- u2 <- 0.7
    b2 <- c(-0.1, -0.1, 0.3, 0.3, 0.6, 0.1, -0.1, -0.8)
    
    w_all <- w_all_2 <- c(-0.1, -0.1, w, -0.1)
    
  }else if ( spp == 'PSSP'){
    a_mu <- a_mu2 <- 2.7
    sig_a <- sig_a2 <- 1.3
    sig_G <- sig_G2 <-0.4
    w <- w2 <-  -1.8
    a <- rep(a_mu, nyrs)
    a2 <- rep(a_mu, nyrs2)
    gint <- gint2 <- rep( 0, G)
    theta <- theta2 <- 1.23
    u <- u2 <- 0.8
    b2 <- c(-0.1, -0.6, 0.5, 0.4, 0.4, 0.4, 0.1, -0.5)
    
    w_all <- w_all_2 <-c(0.16, -0.35, -0.04, w)
    
  }

  inits <- list( NA ) 
  inits[[1]] <- list(a_mu = a_mu, sig_a = sig_a, sig_G = sig_G, a = a, gint = gint, theta = theta, u = u , w = w, 
                     a_mu2 = a_mu2, sig_a2 = sig_a2, sig_G2 = sig_G2, a2 = a2, gint2 = gint2, theta2 = theta2, u2 = u2 , w2 = w2) 
  inits[[2]] <- list(a_mu = a_mu, sig_a = sig_a, sig_G = sig_G, a = a, gint = gint, theta = theta, u = u , w = w, b2 = b2, 
                     a_mu2 = a_mu2, sig_a2 = sig_a2, sig_G2 = sig_G2, a2 = a2, gint2 = gint2, theta2 = theta2, u2 = u2 , w2 = w2)
  inits[[3]] <- list(a_mu = a_mu, sig_a = sig_a, sig_G = sig_G, a = a, gint = gint, theta = theta, u = u , w = w_all, b2 = b2, 
                     a_mu2 = a_mu2, sig_a2 = sig_a2, sig_G2 = sig_G2, a2 = a2, gint2 = gint2, theta2 = theta2, u2 = u2 , w2 = w_all_2)
  
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


