# ##################################################################################################################
# #
# # Get init vals for stan models
# #
# ##################################################################################################################
# 
# #---model descriptions --------------------------------------------------------------------------------------------
# #  
# #   m1: null model:                       w/ size  + intra-specific competition                      
# #   m2: single species model              w/ size  + intra-specific competition + climate 
# #   m3: multi species model               w/ size  + all competition + climate  
# #
# # ----------------------------------------------------------------------------------------------------------------- 
# 
# rm(list = ls())
# 
# library(lme4)
# library(parallel)
# 
# make_df <- function( x ) { 
#   N <- x$N
#   lens <- lapply( x, length )
#   nrs <- lapply(  x, nrow ) 
#   nrs [ sapply( nrs, is.null )  ]  <- 0
#   
#   data.frame( x [ which(lens == N | nrs == N) ]  )
# } 
# 
# set_init_vals_list <-  function( model, C_names, W_names, nyrs2) {  
#   
#   Nspp <- length(W_names[-grep(':', W_names)])
#   init_vals <- as.list( fixef(model)[1:2] )
#   
#   init_vals <- c(init_vals, as.numeric(data.frame( VarCorr(model) )$sdcor[ c(1,2,4)]))
#   
#   names( init_vals )[1:5] <- c('a_mu', 'b1_mu', 'sig_a', 'sig_b1', 'sig_G')
#   
#   if(init_vals$sig_G < 0.0001) { init_vals$sig_G <- 0.05 } # prevent sig G from being zero in some cases 
#   
#   b2 <- fixef(model)[ names(fixef(model))[grep('^[PT]', names(fixef(model)))] ]
#   w <- fixef(model)[ names(fixef(model))[grep('^[W]', names(fixef(model))) ] ]
# 
#   if( length(b2[!is.na(b2)]) > 0 )  { 
#     init_vals <- c(init_vals, list( b2 = as.numeric(b2[!is.na(b2)])))
#   } 
#   if( length(w[!is.na(w)]) > 0 )  { 
#     init_vals <- c(init_vals, list( w = as.numeric(w[!is.na(w)])))  
#   }
#   
#   # random effects initial values need to be specified for STAN 2.6.0
#   nyrs <- length(unique(model@frame$yid) )
#   G <- length(unique(model@frame$gid)) 
#   
#   init_vals$a <- rnorm(nyrs, 0, 0.001)
#   init_vals$b1 <- rnorm(nyrs, 0, 0.001) 
#   init_vals$gint <- rnorm(G, 0, 0.001) 
#   
#   # initial values for the year effects model 
#   init_vals$a2 <- rnorm(nyrs2, 0, 0.001)
#   init_vals$b12 <- rnorm(nyrs2, 0, 0.001)
#   init_vals$gint2 <- init_vals$gint
#   init_vals$w2 <- init_vals$w
#   init_vals$a_mu2 <- init_vals$a_mu
#   init_vals$b1_mu2 <- init_vals$b1_mu
#   init_vals$sig_a2 <- init_vals$sig_a
#   init_vals$sig_b12 <- init_vals$sig_b1
#   init_vals$sig_G2 <- init_vals$sig_G
#   
#   return( init_vals )
# 
# }
# 
# 
# get_init_vals <- function( spp, df, ... ) {
#   
#   nyrs2 <- length(unique(df$treat_year))
#   
#   df <- subset(df, Period == 'Historical')
#   
#   C_names <- names(df)[ grep('^[TP]\\.', names(df))] # climate effects 
#   W_names <- names(df)[ grep('^W', names(df))] # competition effects 
#   W_intra <- names(df)[ grep(spp, names(df))]
#   
#   # write models --------------------------------------------------------------- # 
#   f0 <- 'Y ~ X + (1|gid) + (X|yid)'
#   f1 <- paste(f0, paste(W_intra, collapse = ' + ' ), sep = ' + ')
#   f2 <- paste(f0, paste(c(W_intra, C_names), collapse = ' + '), sep = ' + ')
#   f3 <- paste(f0, paste(c(W_names, C_names), collapse = ' + '), sep = ' + ')
#   f4 <- paste(f0, paste(W_names, collapse = ' + '), sep = ' + ')
#   
#   fs <- list(f1, f2, f3)#, f4) 
#   # ----------------------------------------------------------------------------- #
#   
#   ms <- mclapply( fs, FUN = function( x, ... ) glmer( x , data = df, family = 'binomial'), mc.cores = 4 ) # run models 
#   
#   # set initial values ----------------------------------------
#   init_vals <- lapply( ms, set_init_vals_list, C_names = C_names, W_names = W_names , nyrs2 = nyrs2) 
#   
#   return(init_vals)
# }
# 
# # input files ----------------------------------------------------------------------#
# 
# dfs <- lapply( dir( 'data/temp_data/', '*scaled_survival_dataframe.RDS', full.names = T), readRDS)
# 
# dfs <- dfs[1]
# 
# #nchains <- 4
# spp <- unlist( lapply( dfs, function(x) unique(x$species)) ) 
# spp <- 'ARTR'
# # run functions---------------------------------------------------------------------# 
# init_vals <- mapply( get_init_vals, spp = spp , df = dfs, USE.NAMES = TRUE, SIMPLIFY = FALSE)
# 
# # save output ----------------------------------------------------------------------#
# 
# names(init_vals) <- spp 
# 
# saveRDS( init_vals, file = file.path('data', 'temp_data', 'survival_init_vals.RDS'))
# 
# #

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
#   
# ----------------------------------------------------------------------------------------------------------------- 

rm(list = ls())

get_init_vals_recruitment_models <- function( spp, df, ... ) {
  
  df_old <- subset(df, Period == 'Historical')
  
  nyrs2 <- length(unique(df$treat_year))
  
  nyrs <- length(unique(df_old$yid))
  G <- length(unique(df$gid))
  Ccovs <- length(df[, grep('^[PT]\\.', names(df))])
  
  if(spp == 'ARTR'){
    a_mu <- a_mu2 <- -2.4
    sig_a <- sig_a2 <- 0.6
    sig_G <- sig_G2 <- 0.7
    b1_mu <- 1.7
    sig_b1 <- 0.16
    w <- w2 <- -1.6
    a <- rep(a_mu, nyrs)
    a2 <- rep(a_mu, nyrs2)
    b1 <- rep(b1_mu, nyrs)
    b12 <- rep(b1_mu, nyrs2)
    gint <- gint2 <- rep( 0, G)
    b2 <- rep(0, Ccovs)
    w_all <- w_all_2 <- c(w, -0.1, -0.1, 0, 0, 0.18)
    
  }else if(spp == 'HECO'){
    a_mu <- a_mu2 <- 2.4
    sig_a <- sig_a2 <- 0.6
    sig_G <- sig_G2 <- 0.7
    b1_mu <- 1.7
    sig_b1 <- 0.16
    w <- w2 <- -1.6
    a <- rep(a_mu, nyrs)
    a2 <- rep(a_mu, nyrs2)
    b1 <- rep(b1_mu, nyrs)
    b12 <- rep(b1_mu, nyrs2)
    gint <- gint2 <- rep( 0, G)
    b2 <- c(-0.2, 0.1, 0.2, -0.18, -0.2, 0.2, -0.2, -1.8, -0.4, -0.8, 0.16, 0.15, 0.09, -0.15)
    w_all <- w_all_2 <- c(-0.1, w, -0.1, 0, 0, 0.18)
    
  }else if ( spp == 'POSE'){
    a_mu <- a_mu2 <- 2.4
    sig_a <- sig_a2 <- 0.6
    sig_G <- sig_G2 <- 0.7
    b1_mu <- 1.7
    sig_b1 <- 0.16
    w <- w2 <- -1.6
    a <- rep(a_mu, nyrs)
    a2 <- rep(a_mu, nyrs2)
    b1 <- rep(b1_mu, nyrs)
    b12 <- rep(b1_mu, nyrs2)
    gint <- gint2 <- rep( 0, G)
    b2 <- c(-0.2, 0.1, 0.2, -0.18, -0.2, 0.2, -0.2, -1.8, -0.4, -0.8, 0.16, 0.15, 0.09, -0.15)
    w_all <- w_all_2 <- c(-0.1, -0.1, w, 0, 0, 0.18)
    
  }else if ( spp == 'PSSP'){
    a_mu <- a_mu2 <- 2.4
    sig_a <- sig_a2 <- 0.6
    sig_G <- sig_G2 <- 0.7
    b1_mu <- 1.7
    sig_b1 <- 0.16
    w <- w2 <- -1.6
    a <- rep(a_mu, nyrs)
    a2 <- rep(a_mu, nyrs2)
    b1 <- rep(b1_mu, nyrs)
    b12 <- rep(b1_mu, nyrs2)
    gint <- gint2 <- rep( 0, G)
    b2 <- c(-0.2, 0.1, 0.2, -0.18, -0.2, 0.2, -0.2, -1.8, -0.4, -0.8, 0.16, 0.15, 0.09, -0.15)
    w_all <- w_all_2 <- c(-0.1, -0.1, 0, w, 0, 0.18)
    
  }
  
  inits <- list( NA ) 
  inits[[1]] <- list(a_mu = a_mu, b1_mu = b1_mu, sig_b1 = sig_b1,  sig_a = sig_a, sig_G = sig_G, a = a, gint = gint, w = w, 
                     a_mu2 = a_mu2, b1_mu = b1_mu, sig_b1 = sig_b1, sig_a2 = sig_a2, sig_G2 = sig_G2, a2 = a2, gint2 = gint2, w2 = w2) 
  inits[[2]] <- list(a_mu = a_mu, b1_mu = b1_mu, sig_b1 = sig_b1, sig_a = sig_a, sig_G = sig_G, a = a, gint = gint, w = w, b2 = b2, 
                     a_mu2 = a_mu2, b1_mu = b1_mu, sig_b1 = sig_b1, sig_a2 = sig_a2, sig_G2 = sig_G2, a2 = a2, gint2 = gint2, w2 = w2)
  inits[[3]] <- list(a_mu = a_mu, b1_mu = b1_mu, sig_b1 = sig_b1, sig_a = sig_a, sig_G = sig_G, a = a, gint = gint, w = w_all, b2 = b2, 
                     a_mu2 = a_mu2, b1_mu = b1_mu, sig_b1 = sig_b1, sig_a2 = sig_a2, sig_G2 = sig_G2, a2 = a2, gint2 = gint2, w2 = w_all_2)
  
  return(inits)
}

# input files ----------------------------------------------------------------------#

dfs <- lapply( dir( 'data/temp_data/', '*scaled_survival_dataframe.RDS', full.names = T), readRDS)
spp <- unlist( lapply( dfs, function(x) unique(x$species)) ) 

dfs <- dfs[1]
spp <- 'ARTR'
# run functions---------------------------------------------------------------------# 
init_vals <- list( NA )

init_vals <- mapply(get_init_vals_recruitment_models, spp = factor(spp), df = dfs, SIMPLIFY = FALSE)

names(init_vals) <- spp 

# save output ----------------------------------------------------------------------#

saveRDS( init_vals, file = file.path('data/temp_data', 'survival_init_vals.RDS'))


