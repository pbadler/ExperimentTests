##################################################################################################################
#
# Fit models to historical data only 
#
# Export starting values
#
##################################################################################################################

#---model descriptions --------------------------------------------------------------------------------------------
#  
#   m0: null model:                       
#   m1: climate model:                    w/ climate 
#   m2: single species model:             w/ climate, w/ intra-specific competition 
#   m3: single species climate model:     w/ climate, w/ intra-specific competition
#   m4: full model,                       w/ climate, w/ intra-specific competition, w/ interspecific competition 
#
# ----------------------------------------------------------------------------------------------------------------- 

rm(list = ls())

library(lme4)

dl <- readRDS('data/temp_data/data_lists_for_stan_models.RDS')

test <- dl[[1]]
spp <- names(dl)[1]
nchains <- 4 

make_df <- function( sl ) { 
  
  N <- sl$N
  lens <- lapply( sl, length )
  nrs <- lapply(  sl, nrow ) 
  nrs [ sapply( nrs, is.null )  ]  <- 0
  
  data.frame( sl [ which(lens == N | nrs == N) ]  )
} 

test <-  make_df(test ) 


Cnames <- names(test)[ grep('^C', names(test))] # climate effects 
Wnames <- names(test)[ grep('^W', names(test))] # competition effects 
Wintra <- names(test)[ grep(spp, names(test))]

# write models --------------------------------------------------------------- # 
f0 <- 'Y ~ X + (1|gid) + (X|yid)'
f1 <- paste(f0, paste(Cnames, collapse = ' + ' ), sep = ' + ')
f2 <- paste(f0, paste(Wintra, collapse = ' + ' ), sep = ' + ')
f3 <- paste(f0, paste(c(Wintra, Cnames), collapse = ' + '), sep = ' + ')
f4 <- paste(f0, paste(c(Wnames, Cnames), collapse = ' + '), sep = ' + ')

fs <- c(f0, f1, f2, f3, f4 ) 
# ----------------------------------------------------------------------------- #

ms <- lapply(fs, FUN = function( x, ... ) lmer( x, data = test) ) # run models 

m_null <- ms[[5]]

refx <- as.data.frame(VarCorr(m_null))$sdcor
init_vals <- c(as.numeric(fixef(m_null)), refx)

init_vals <- cbind(  c(names( fixef(m_null)), as.data.frame(VarCorr(m_null))$grp), init_vals  )

init_vals

# set initial values ----------------------------------------

b2 <- fixef(m_null)[Cnames]
w <- fixef(m_null)[Wnames]



inits <- rep ( list(list( sigma = init_vals[nrow(init_vals), 2 ],
                          a_mu = init_vals[1, 2],
                          b1_mu = init_vals[ 2, 2 ]  ) ), nchains  )

inits
