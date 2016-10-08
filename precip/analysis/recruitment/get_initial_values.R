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
  
  C_names <- names(df)[ grep('^[PT]\\.', names(df))] # climate effects 
  parents1_names <- names(df)[ grep('^cov', names(df))] # competition effects 
  parents2_names <- names(df)[ grep('^Gcov', names(df))] # competition effects 
  
  parents_intra <- names(df)[ grep(spp[1], names(df))]
  
  u <- 0.9 ### mixing parameter 
  
  # need to do some transformation to get parameters correct ------------------------ # 
  
  effective_cover <- df[, parents1_names ]*u + df[, parents2_names]*(1-u)
  
  effective_cover_sqrd <- sqrt( effective_cover )
  
  intra_df <- data.frame( 'intra' = log( effective_cover[, parents_intra[1]])) 
  
  prepped_df <- data.frame( df[ , c(C_names, 'Y', 'gid', 'yid')], effective_cover_sqrd, intra_df )
  
  W_names <- names(prepped_df)[grep('^cov\\.', names(prepped_df))]
  W_intra <- names(prepped_df)[grep(spp, names(prepped_df))]
  
  prepped_df$Y <- df$Y
  
  # write models --------------------------------------------------------------- # 
  f0 <- paste( 'Y', ' ~ offset(intra)', '+ (1|gid) + (1|yid)' )
  f1 <- paste(f0, paste(W_intra, collapse = ' + ' ), sep = ' + ')
  f2 <- paste(f0, paste(c(W_intra, C_names), collapse = ' + '), sep = ' + ')
  f3 <- paste(f0, paste(c(W_names, C_names), collapse = ' + '), sep = ' + ')
  f4 <- paste(f0, paste(W_names, collapse = ' + '), sep = ' + ')
  
  fs <- list(f1, f2, f3, f4 )
  # ----------------------------------------------------------------------------- #
  
  ms <- mclapply( fs, FUN = function( x, ... ) glmer.nb( x , data = prepped_df) , mc.cores = 4) # run models 
  
  # set initial values ----------------------------------------
  init_vals <- lapply( ms, set_init_vals_list, C_names = C_names, W_names = W_names ) 
   
  return(init_vals)
}

# input files ----------------------------------------------------------------------#
 
dfs <- lapply( dir( 'data/temp_data/', '*scaled_recruitment_dataframe.RDS', full.names = T), readRDS)
spp <- unlist( lapply( dfs, function(x) unique(x$species)) ) 

# run functions---------------------------------------------------------------------# 
init_vals <- list( NA )

init_vals[[1]] <- get_init_vals_recruitment_models( spp[1], dfs[[1]])
init_vals[[2]] <- get_init_vals_recruitment_models( spp[2], dfs[[2]]) ### Throwing error 
init_vals[[3]] <- get_init_vals_recruitment_models( spp[3], dfs[[3]])
init_vals[[4]] <- get_init_vals_recruitment_models( spp[4], dfs[[4]])

names(init_vals) <- spp 
# save output ----------------------------------------------------------------------#

saveRDS( init_vals, file = file.path('data/temp_data', 'recruitment_init_vals.RDS'))


