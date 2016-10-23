rm(list = ls())

load('data/temp_data/master_list.Rdata')

df <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')[['HECO']]
sdf <- readRDS('data/temp_data/survival_data_lists_for_stan.RDS')[['HECO']]

library(lme4)
spp_df <- readRDS('data/temp_data/HECO_scaled_growth_dataframe.RDS')
spp_df <- subset(spp_df, Period == 'Historical')
m1 <- lmer( data = spp_df, Y ~ gid + X + W.HECO + (X|yid)) 
summary(m1)
ranef(m1)

b1_mu = fixef(m1)['X']
b1 = ranef(m1)[[1]][, 2] + b1_mu

inits <- list( a = ranef(m1)[[1]][, 1], b1_mu = b1_mu, b1 = b1, bg = fixef(m1)[1:df$G], w = fixef(m1)['W.HECO'], sig_a = 1.5, sig_b1 = 0.25, tauSize = -0.5, tau = 4)

inits2 <- readRDS('data/temp_data/growth_init_vals.RDS')[['HECO']]

library(rstan)
df$tau_beta <- as.numeric( master_list$sd_vec[2, 2] ) 
df$W <- df$W[, df$spp]
df$W2 <- df$W2[, df$spp]
df$Whold <- df$Whold[, df$spp]

testfit <- stan('analysis/growth/model_growth_1.stan', chains = 4, data = df, init = rep(list(inits2[[1]]), 4) , cores = 4, iter = 100)

traceplot(testfit, 'w')
traceplot( testfit, 'bg')
traceplot( testfit, 'b1_mu')
traceplot(testfit, 'sig_a')
traceplot(testfit, 'tau')
traceplot(testfit, 'tauSize')
summary(testfit, c('sig_a'))$summary
summary(testfit, c('b1_mu', 'tau', 'tauSize', 'bg', 'sig_a', 'sig_b1', 'w'))$summary[, 1]

plot(ranef(m1)[[1]][, 1], type ='l')
points( summary(testfit, c('a'))$summary[, 6], type = 'l', col = 'red')

plot(ranef(m1)[[1]][, 2], type ='l')
points( summary(testfit, c('b1'))$summary[, 6] - summary(testfit, 'b1_mu')$summary[, 1], type = 'l', col = 'red') 

test_preds <- rstan::extract(testfit, c('a', 'b1'))

test_df <- data.frame( test_preds ) 
library(dplyr)
library(tidyr)
library(stringr)

b1_df<- test_df %>% gather( var, val , starts_with('b1')) %>% mutate( var_label = as.numeric(str_extract(pattern = '[0-9]+$', var)), m_b1_mu2 = mean(b1_mu2))

a_df<- test_df %>% gather( var, val , starts_with('a2')) %>% mutate( var_label = as.numeric(str_extract(pattern = '[0-9]+$', var)))

ggplot (b1_df, aes( x = val)) + geom_histogram() + geom_vline(aes( xintercept = m_b1_mu2)) + facet_wrap( ~ var_label, ncol = 4)

ggplot (a_df, aes( x = val)) + geom_histogram() + geom_vline(aes( xintercept = 0)) + facet_wrap( ~ var_label, ncol = 4)


testframe <- readRDS('data/temp_data/HECO_scaled_survival_dataframe.RDS')
testframe$treat_year_label
testframe %>% distinct( yid , treat_year_label) %>% arrange( yid )


hist( test_preds$mu - df$Y )
hist( test_preds$muhat - df$Yhold )
hist( test_preds$muhat2 - df$Yhold)

source( 'analysis/waic_fxns.R')

waic(testfit, llname = 'log_lik2')

summary(testfit, c('bg', 'b1_mu', 'sig_a', 'sig_b1'))$summary

test_df <- data.frame(gid = factor(df$gid), Y = df$Y, yid = df$yid, X = df$X, W.HECO = df$W[,1], P.f.w.sp.l = df$C[, 3], `P.f.w.sp.l:logarea` = df$C[,10])
df$gid

mymod <- glmer( data = test_df, Y ~ X + gid + (X| yid ), family = 'binomial')
summary(mymod)

hist(rstan::extract(testfit,'sig_a')$sig_a)
hist(rstan::extract(testfit, 'sig_b1')$sig_b1)
hist(rstan::extract(testfit, 'b2[14]')$`b2[14]`)

traceplot(testfit, 'bg')
traceplot(testfit, 'b2')
traceplot(testfit, 'b1_mu')
traceplot(testfit, 'sig_b1')
traceplot(testfit, 'sig_a')
traceplot(testfit, 'a')

colMeans(df$C)