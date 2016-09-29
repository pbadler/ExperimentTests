##################################################################################################################
#
# Get init vals for stan models
#
##################################################################################################################

#---model descriptions --------------------------------------------------------------------------------------------
#  
#   m1: null model:                       
#   m2: climate model:                    w/ climate 
#   m3: single species model:                         w/ intra-specific competition 
#   m4: single species climate model:     w/ climate, w/ intra-specific competition
#   m5: full model:                       w/ climate, w/ intra-specific competition, w/ interspecific competition 
#
# ----------------------------------------------------------------------------------------------------------------- 

rm(list = ls())

library(lme4)

make_df <- function( x ) { 
  N <- x$N
  lens <- lapply( x, length )
  nrs <- lapply(  x, nrow ) 
  nrs [ sapply( nrs, is.null )  ]  <- 0
  
  data.frame( x [ which(lens == N | nrs == N) ]  )
} 

set_init_vals_list <-  function( model, C_names, W_names, dl ) {  
  
  init_vals <- as.list( fixef(model)[1:2] )
  
  init_vals <- c(init_vals, as.numeric(data.frame( VarCorr(model) )$sdcor[ c(1,2,4)]))
  
  names( init_vals )[1:5] <- c('a_mu', 'b1_mu', 'sig_a', 'sig_b1', 'sig_G')
  
  if(init_vals$sig_G < 0.0001) { init_vals$sig_G <- 0.05 } # prevent sig G from being zero in some cases 
  
  b2 <- fixef(model)[C_names]
  w <- fixef(model)[W_names]

  if( length(b2[!is.na(b2)]) > 0 )  { 
    init_vals <- c(init_vals, list( b2 = as.numeric(b2[!is.na(b2)])))
  } 
  if( length(w[!is.na(w)]) > 0 )  { 
    init_vals <- c(init_vals, list( w = as.numeric(w[!is.na(w)])))  
  }
  
  # random effects initial values need to be specified for STAN 2.6.0
  nyrs <- length(unique(model@frame$yid) )
  G <- length(unique(model@frame$gid)) 
  
  init_vals$a <- rep(0, nyrs)
  init_vals$b1 <- rep(0, nyrs) 
  init_vals$gint <- rep(0, G) 
  
  return( init_vals )

}


get_init_vals_survival_models <- function( spp, datalist, ... ) {

  df <- make_df(x = datalist) 
  
  C_names <- names(df)[ grep('^C', names(df))] # climate effects 
  W_names <- names(df)[ grep('^W', names(df))] # competition effects 
  W_intra <- names(df)[ grep(spp, names(df))]
  
  # write models --------------------------------------------------------------- # 
  f1 <- 'Y ~ X + (1|gid) + (X|yid)'
  f2 <- paste(f1, paste(C_names, collapse = ' + ' ), sep = ' + ')
  f3 <- paste(f1, paste(W_intra, collapse = ' + ' ), sep = ' + ')
  f4 <- paste(f1, paste(c(W_intra, C_names), collapse = ' + '), sep = ' + ')
  f5 <- paste(f1, paste(c(W_names, C_names), collapse = ' + '), sep = ' + ')
  
  fs <- list(f1, f2, f3, f4, f5 ) 
  # ----------------------------------------------------------------------------- #
  
  ms <- lapply( fs, FUN = function( x, ... ) glmer( x , data = df, family = 'binomial') ) # run models 
  
  # set initial values ----------------------------------------
  init_vals <- lapply( ms, set_init_vals_list, C_names = C_names, W_names = W_names ) 
  
  return(init_vals)
}

# input files ----------------------------------------------------------------------#

dl <- readRDS('data/temp_data/survival_data_lists_for_stan.RDS')
nchains <- 4
spp <- names(dl)

# run functions---------------------------------------------------------------------# 
init_vals <- mapply( get_init_vals_survival_models, spp = spp , datalist = dl, USE.NAMES = TRUE, SIMPLIFY = FALSE)

# save output ----------------------------------------------------------------------#

saveRDS( init_vals, file = file.path('data', 'temp_data', 'survival_init_vals.RDS'))


