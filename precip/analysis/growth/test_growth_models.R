rm(list  = ls())
library(rstan)
library(lme4)
# simulate climate, competition, year and group effects ------------------------------------- # 

test_dat <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[['ARTR']]

sig_a <- 1
sig_b1 <- 0.5
b1_mu  <- 1.4

pars <- list( 
  bg = c(1.5, runif(test_dat$G - 1, -1, 1)),
  b2 = c(2, rep(0, test_dat$Covs - 1)),
  a  = rnorm(test_dat$nyrs, 0, sd = sig_a),
  b1 = rnorm(test_dat$nyrs, b1_mu, sd = sig_b1),
  w  = c(1, 0, 0, 0),
  sigma = 1)

simulate_growth <- function( pars , test_dat ){ 
  
  N <- test_dat$N
  C <- test_dat$C
  W <- test_dat$W
  X <- test_dat$X
  gm <- test_dat$gm
  yid <- test_dat$yid

  attach(pars )
  
  climEff  <- C%*%b2
  gint     <- gm%*%bg
  coverEff <- W%*%w
  
  mu <- coverEff
  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + coverEff[n] + climEff[n]
  }

  Y <- rnorm(N, mu, sigma)
  rm(pars)
  detach(pars )   
  return(Y)
}

test_dat$Y <- simulate_growth(pars, test_dat)

myfit1 <- stan('analysis/growth/model_growth_1.stan', data = test_dat, chains = 4, cores = 4, iter = 2000)

traceplot( myfit1, 'b1')
traceplot(myfit1, 'b1_mu')
traceplot(myfit1, 'sig_b1')
traceplot( myfit1, 'b2')
traceplot( myfit1, 'w')

my_fit_summary <- summary(myfit1, c('a', 'sig_a', 'b1', 'b1_mu', 'sig_b1', 'sigma', 'b2', 'w', 'bg'))$summary

my_fit_summary <- data.frame( my_fit_summary ) 
my_fit_summary$par <- row.names( my_fit_summary)

df <- cbind( par = pars$a, est = my_fit_summary[1:21, 1] )
plot(df)

plot( cbind( par = pars$b1, est = my_fit_summary[23:43, 1] ) )
abline(0,1)

mydf <- as.data.frame( test_dat[c('X' , 'Y', 'yid', 'gid')] )
mydf$W <- test_dat$W
mydf$C <- test_dat$C
m1 <- lmer( data = mydf, Y ~ X + factor(gid) + (X|yid) + C + W  )

plot( cbind( my_fit_summary[23:43, 1] - b1_mu,  ranef(m1)$yid[,2]) )

plot(pars$a, ranef(m1)$yid[,1])
plot(pars$b1 - b1_mu, ranef(m1)$yid[,2])
abline(0,1)

# simulate year effects 


sig_a <- 0.5
sig_b1 <- 0.1
b1_mu  <- 0.8

pars <- list( 
  bg = c(1.5, runif(test_dat$G - 1, -1, 1)),
  b2 = c(0, rep(0, test_dat$Covs - 1)),       # no climate effects 
  a  = rnorm(test_dat$nyrs, 0, sd = sig_a),
  b1 = rnorm(test_dat$nyrs, b1_mu, sd = sig_b1),
  w  = c(1, 0, 0, 0),
  sigma = 1)

test_dat$Y <- simulate_growth(pars, test_dat)

myfit1 <- stan('analysis/growth/model_growth_1.stan', data = test_dat, chains = 4, cores = 4, iter = 2000)

my_fit_summary <- summary(myfit1, c('a', 'sig_a', 'b1', 'b1_mu', 'sig_b1', 'sigma', 'b2', 'w', 'bg'))$summary

my_fit_summary <- data.frame( my_fit_summary ) 
my_fit_summary$par <- row.names( my_fit_summary)

df <- cbind( par = pars$a, est = my_fit_summary[1:21, 1] )
plot(df)

df <- cbind( par = pars$b1, est = my_fit_summary[ grep( 'b1\\_' , my_fit_summary$par), 1 ] )



# simulate treatment effects 
rm(list = ls())
test_dat <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[['ARTR']]

sig_a <- 1
sig_b1 <- 0.6
b1_mu  <- 1.5

pars <- list( 
  bg = c(1.5, runif(test_dat$G - 1, -1, 1)),
  bt = c(1, -1, 1, -1),
  a  = rnorm(test_dat$nyrshold, 0, sd = sig_a),
  b1 = rnorm(test_dat$nyrshold, b1_mu, sd = sig_b1),
  w  = c(1, 0, 0, 0),
  sigma = 0.5)

simulate_growth_treatment <- function( pars , test_dat ){ 
  nyrs <- test_dat$nyrs
  N <- test_dat$Nhold
  C <- test_dat$Chold
  W <- test_dat$Whold
  X <- test_dat$Xhold
  gm <- test_dat$gmhold
  yid <- test_dat$yidhold
  tm <- test_dat$tmhold
  
  attach(pars )
  
  treatEff  <- tm%*%bt
  gint     <- gm%*%bg
  coverEff <- W%*%w
  
  mu <- coverEff
  
  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n] -nyrs] + b1[yid[n] - nyrs]*X[n] + coverEff[n] + treatEff[n]
  }
  
  Y <- rnorm(N, mu, sigma)
  rm(pars)
  detach(pars )   
  return(Y)
}

# treatment model 
test_dat$Yhold <- simulate_growth_treatment(pars, test_dat)

myfit1 <- stan('analysis/growth/model_growth_treatment_effects.stan', data = test_dat, chains = 4, cores = 4, iter = 2000 )

a <- summary(myfit1, c('a'))$summary[, 1]
b1 <- summary(myfit1, c('b1'))$summary[,1]
b1_mu_est <- summary(myfit1, c('b1_mu'))$summary[,1]
b1_mu 
b1_mu_est

plot ( pars$a, a )
abline(0,1)

plot( pars$b1, b1 )
abline(0,1)

bt <- summary( myfit1 , c('bt'))$summary
bt
pars$bt

mydf <- as.data.frame( test_dat[c('Xhold' , 'Yhold', 'yidhold', 'gidhold', 'treathold')] )
mydf$W <- test_dat$Whold
m1 <- lmer( data = mydf, Yhold ~ Xhold + W + factor(treathold) + factor(gidhold) + (Xhold|yidhold)  )
summary(m1)

plot( pars$a, ranef(m1)$yidhold[,1] )
abline(0,1)
plot( pars$b1 , ranef(m1)$yidhold[,2] + b1_mu )
abline(0, 1)

plot( a, ranef(m1)$yidhold[,1])
abline(0,1)
plot( b1, ranef(m1)$yidhold[,2] + b1_mu)
abline(0,1)
