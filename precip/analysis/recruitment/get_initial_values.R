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

set_init_vals_list <-  function( model, C_names, W_names ) {  
  
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
  
  nyrs <- length(unique(df$yid))
  G <- length(unique(df$gid))
  Ccovs <- length(df[, grep('^[PT]\\.', names(df))])
  
  a_mu <- 4.0
  sig_a <- 1.5
  sig_G <- 0.2
  w <- -3.4
  a <- rep(a_mu, nyrs)
  gint <- rep( 0, G)
  theta <- 0.6
  u <- 0.25
  b2 <- rep(0, Ccovs)
  
  w2 <- w*diag(4)
  w2[w2==0] <- 0.01
  
  inits <- list( NA ) 
  
  inits[[1]] <- list(a_mu = a_mu, sig_a = sig_a, sig_G = sig_G, a = a, gint = gint, theta = theta, u = u , w = w)
  inits[[2]] <- list(a_mu = a_mu, sig_a = sig_a, sig_G = sig_G, a = a, gint = gint, theta = theta, u = u , w = w, b2 = b2)
  inits[[3]] <- list(a_mu = a_mu, sig_a = sig_a, sig_G = sig_G, a = a, gint = gint, theta = theta, u = u , w = w2[spp, ], b2 = b2)

  return(inits)
}

# input files ----------------------------------------------------------------------#
 
dfs <- lapply( dir( 'data/temp_data/', '*scaled_recruitment_dataframe.RDS', full.names = T), readRDS)
spp <- unlist( lapply( dfs, function(x) unique(x$species)) ) 

# run functions---------------------------------------------------------------------# 
init_vals <- list( NA )

dfs <- lapply( dfs, function(x) subset( x , Period == 'Historical'))

init_vals <- mapply(get_init_vals_recruitment_models, spp = factor(spp), df = dfs, SIMPLIFY = FALSE)

names(init_vals) <- spp 

# save output ----------------------------------------------------------------------#

saveRDS( init_vals, file = file.path('data/temp_data', 'recruitment_init_vals.RDS'))


