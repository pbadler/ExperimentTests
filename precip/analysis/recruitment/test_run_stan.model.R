rm(list = ls())
library(rstan)

datalist <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')

test_spp <- 'POSE'
test <- datalist[[test_spp]]

test$parents1 <- test$parents1 [ , grep( colnames (test$parents1), pattern = test_spp)  ] 
test$parents2 <- test$parents2 [ , grep( colnames (test$parents2), pattern = test_spp)  ] 

test$parents2[test$parents2 == 0 ] <- min(test$parents2[test$parents2 > 0 ])

inits <- readRDS('data/temp_data/recruit_init_vals.RDS')


inits <- inits[[3]]
inits
inits <- rep(list( inits), 3)

test$tau_beta <- 1.5

test <- test[c('N', 'Y', 'gid', 'G', 'Yrs', 'yid', 'C', 'Covs', 'parents1', 'parents2', 'tau_beta')]

test_fit <- stan('analysis/recruitment/simple_recruitment.stan', 
                 data = test, init = inits, chains = 3, iter = 2000)

test_fit

traceplot(test_fit, c('dd', 'b2'))
traceplot(test_fit, c('theta', 'u'))
traceplot(test_fit, c('a_mu'))

# test model 2 ------------------------------------------------------------------------------------
rm(list = ls())

library(rstan)

datalist <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')

test_spp <- 'POSE'
test <- datalist[[test_spp]]
test$spp

test$parents2[test$parents2 == 0 ] <- min(test$parents2[test$parents2 > 0 ])

inits <- list(a = rep( 3, test$Yrs),
              a_mu = 3.4, 
              sig_a = 0.6,
              gint = rep( 0, test$G), 
              sig_G = 0.15, 
              u = 0.5, 
              dd = rep(-1, test$Nspp), 
              theta=1.3)
inits <- rep(list(inits), 1)

test$tau_beta <- 1.5

test <- test[c('N', 'Nspp', 'spp',  'Y', 'gid', 'G', 'Yrs', 'yid', 'C', 'Covs', 'parents1', 'parents2', 'tau_beta')]

test_fit <- stan('analysis/recruitment/simple_recruitment2.stan', 
                 data = test, init = inits, chains = 1, iter = 1000, cores = 1)


traceplot(test_fit, c('b2'))
traceplot(test_fit, c('dd'))
traceplot(test_fit, c('theta'))
traceplot( test_fit, c('a_mu'))
traceplot( test_fit, 'u')


# test with predictions --------------------------------------------------------------------------------------------------- # 
rm(list = ls())

library(rstan)

datalist <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')

test_spp <- 'POSE'
test <- datalist[[test_spp]]

test$parents2[test$parents2 == 0 ] <- min(test$parents2[test$parents2 > 0 ])

inits <- list(a = rep( 3, test$Yrs),
              a_mu = 3.4, 
              sig_a = 0.6,
              gint = rep( 0, test$G), 
              sig_G = 0.15, 
              u = 0.5, 
              dd = rep(-1, test$Nspp), 
              theta=1.3)

inits <- rep(list(inits), 1)

test$tau_beta <- 1.5

test <- test[c('N', 'Nspp', 'spp',  'Y', 'gid', 'G', 'Yrs', 'yid', 'C', 'Covs', 'parents1', 'parents2', 'tau_beta', 
               'npreds', 'nyrs_out', 'yid_out', 'gid_out', 'y_holdout', 'Chold', 'parents1_out', 'parents2_out')]


test_fit <- stan('analysis/recruitment/simple_recruitment_with_pred.stan', 
                 data = test, init = inits, chains = 1, iter = 100, cores = 1)


traceplot(test_fit, c('b2'))
traceplot(test_fit, c('dd'))
traceplot(test_fit, c('theta'))
traceplot( test_fit, c('a_mu'))
traceplot( test_fit, 'u')

# ----------------------------------------------------------------------------------------------------------------- # 

library(lme4)

test <- datalist$POSE

make_df <- function( x ) { 
  N <- x$N
  lens <- lapply( x, length )
  nrs <- lapply(  x, nrow ) 
  nrs [ sapply( nrs, is.null )  ]  <- 0
  
  data.frame( x [ which(lens == N | nrs == N) ]  )
} 


test_df <- make_df(test)

N <- test$N
lens <- lapply( test, length)
nrs <- lapply( test, nrow)

str(test)


lapply(test2$ARTR, length)
str(test2$ARTR)

test2 <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')

head( make_df(test2$ARTR) ) 

head( test_df ) 

m1 <- glmer.nb(data = test_df, Y ~ C.P.a.0 + (1|gid))

m1@theta
summary(m1)
