rm(list  = ls())
library(rstan)
# simulate climate, competition, year and group effects ------------------------------------- # 

test_dat <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[['POSE']]

sig_a <- 0.4
sig_b1 <- 0.4
b1_mu  <- 0.8

pars <- list( 
  bg = c(1.5, runif(test_dat$G - 1, -1, 1)),
  b2 = c(2, rep(0, test_dat$Covs - 1)),
  a  = rnorm(test_dat$nyrs, 0, sd = sig_a),
  b1 = rnorm(test_dat$nyrs, b1_mu, sd = sig_b1),
  w  = c(1, 0, 0, 0))

simulate_survival <- function( pars , test_dat ){ 
  
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
  
  p <- exp(mu)/(1 + exp(mu)) # inverse logit 
  
  Y <- rbinom(N, size = 1,  p)
  rm(pars)
  detach(pars )   
  return(Y)
}

# test_dat$Y <- simulate_survival(pars, test_dat)
# 
# myfit1 <- stan('analysis/survival/model_survival_1.stan', data = test_dat, chains = 4, cores = 4, iter = 1000, seed = 1)
# 
# ye <- summary(myfit1, c('a', 'b1'))$summary[,1]
# plot( ye[1:21], c(pars$a))
# abline(0,1)
# plot( ye[22:42], pars$b1)
# abline(0,1)
# 
# pars$b2
# summary(myfit1, 'b2')$summary[,1]
#               


# simulate treatment effects 
rm(list = ls())
test_dat <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[['POSE']]

sig_a <- 0.8
sig_b1 <- 0.4
b1_mu  <- 2

pars <- list( 
  bg = c(1.5, runif(test_dat$G - 1, -1, 1)),
  bt = c(1, -1, 1, -1),
  a  = rnorm(test_dat$nyrshold, 0, sd = sig_a),
  b1 = rnorm(test_dat$nyrshold, b1_mu, sd = sig_b1),
  w  = c(1, 0, 0, 0), 
  b1_mu = b1_mu)

simulate_survival_treatment <- function( pars , test_dat ){ 
  nyrs <- test_dat$nyrs
  N <- test_dat$Nhold
  tm <- test_dat$tmhold
  W <- test_dat$Whold
  X <- test_dat$Xhold
  gm <- test_dat$gmhold
  yid <- test_dat$yidhold
  
  attach(pars )
  
  treatEff  <- tm%*%bt
  gint     <- gm%*%bg
  coverEff <- W%*%w
  
  mu <- coverEff
  
  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n] - nyrs] + b1[yid[n] - nyrs]*X[n] + coverEff[n] + treatEff[n]
  }
  
  p <- exp(mu)/(1 + exp(mu)) # inverse logit 
  
  Y <- rbinom(N, size = 1,  p)
  rm(pars)
  detach(pars )   
  return(Y)
}

test_dat$Yhold <- simulate_survival_treatment(pars, test_dat)

myfit1 <- stan('analysis/survival/model_survival_treatment_effects.stan', data = test_dat, chains = 4, cores = 4, iter = 2000 , thin = 4)

b1_mu <- summary(myfit1, 'b1_mu')$summary[,1]
pars$b1_mu

cbind( pars$a, summary(myfit1, 'a')$summary)
cbind( pars$b1, summary(myfit1, 'b1')$summary)

a <- summary(myfit1, 'a')$summary[,1]
plot( pars$a, a)
abline(0,1)

b1 <- summary(myfit1, 'b1')$summary[,1]
plot(pars$b1, b1 )
abline(0,1)
pars$b1

bt <- summary(myfit1, 'bt')$summary[,1]
bt

mu <-summary(myfit1, 'mu')$summary[,1]
plot ( mu, test_dat$Yhold)

w <- summary(myfit1, 'w')$summary[,1]              
w

testdf <- as.data.frame( test_dat[c('yidhold', 'gidhold', 'treathold', 'Xhold', 'Yhold', 'gidhold')])
testdf$W <- test_dat$Whold

library(lme4)
m1 <- glmer( data = testdf, 'Yhold ~ Xhold + W + factor(treathold)*Xhold + (Xhold|yidhold) + factor(gidhold)', family = 'binomial')
summary(m1)

plot( pars$a, ranef(m1)$yidhold[, 1] )
abline(0,1)
plot( pars$b1 - b1_mu, ranef(m1)$yidhold[, 2])
abline(0,1)









source('analysis/waic_fxns.R')

waic(myfit1)

ests1 <- summary(myfit1, c('b1_mu', 'w', 'b2', 'sig_a', 'sig_b1'))$summary[, 1]
ests1

mu <- summary( myfit1, 'mu')$summary[, 1]
muhat <- summary(myfit1, 'muhat')$summary[, 1]

plot( mu, test_dat$Y )
abline(0,1)

plot(muhat, test_dat$Yhold)
abline(0,1)

plot(test_dat$Xhold, muhat)



traceplot( myfit1, 'b2')
traceplot(myfit1, 'w')

# fit same with lmer 
library(lme4)

df <- data.frame( X = test_dat$X, Y = test_dat$Y)
df$C <- test_dat$C
df$W <- test_dat$W
df$Group <- factor(test_dat$gid)
df$yid <- test_dat$yid

f <- as.formula( Y ~ Group + X + (X|yid) + W + C)
m1 <- glmer(data = df, f , family = 'binomial')

summary(m1)

data.frame( summary(m1)$coefficients) 

Ceffects <- summary(myfit1, c('b2'))$summary[, 1]

cbind(fixef(m1)[grep( '^C' , names( fixef(m1)))], Ceffects, real = pars$b2)

# estimates table 

trueVals <- data.frame( type = 'true_val' , b1_mu = b1_mu, sig_b1 = sig_b1, sig_a = sig_a, 
                        Int = pars$bg[1], 
                        w_1 = pars$w[1], w_2 = pars$w[2], w_3 =pars$w[3], w_4 = pars$w[4], 
                        b2_1 = pars$b2[1], b2_2 = pars$b2[2], b2_3 = pars$b2[3], b2_4 = pars$b2[4], 
                        b2_5 = pars$b2[5], b2_6 = pars$b2[6], b2_7 = pars$b2[7], b2_8 = pars$b2[8]) %>% 
  gather(par , set_val , b1_mu:b2_7)

stanVals <- data.frame( summary(myfit1, c('bg[1]', 'b1_mu', 'sig_b1', 'sig_a', 'w', 'b2'))$summary[ , c(1,2)] )
stanVals$par <- row.names(stanVals) 

stanVals$par[ stanVals$par == "bg[1]" ] <- 'Int'
library(stringr)
stanVals$par <- str_replace_all( stanVals$par, pattern = '\\[', replacement = '_')
stanVals$par <- str_replace_all( stanVals$par, pattern = '\\]', replacement = '')

lmer_vals <- data.frame( summary(m1)$coefficients ) 
lmer_vals$par <- row.names(lmer_vals)

lmer_vals$par[ lmer_vals$par == 'X' ]  <- 'b1_mu'
lmer_vals$par[ grep('^WW',  lmer_vals$par) ]  <- paste0( 'w_', c(1:4))
lmer_vals

lmer_vals$par[ grep( '^C', lmer_vals$par ) ] <- paste0('b2_', c(1:8))

lmer_vals$par[lmer_vals$par == '(Intercept)'] <- 'Int'

merge( merge(trueVals, stanVals), lmer_vals, all.x = TRUE)

plot( data = merge( lmer_vals,  stanVals )[-1, ], Estimate ~ mean ) 
