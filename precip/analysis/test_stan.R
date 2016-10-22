rm(list = ls())

load('data/temp_data/master_list.Rdata')

df <- readRDS('data/temp_data/survival_data_lists_for_stan.RDS')[['PSSP']]

library(rstan)

df$tau_beta <- master_list$sd_vec[1, 2]
# df$W <- df$W[, df$spp]
# df$W2 <- df$W2[, df$spp]
# df$Whold <- df$Whold[, df$spp]
testfit <- stan('analysis/survival/model_survival_3.stan', chains = 1, data = df, cores = 1, iter = 1000)

traceplot(testfit, 'b12')

test_preds <- rstan::extract(testfit, c('a', 'a2', 'b1', 'b12', 'b1_mu2'))

test_df <- data.frame( test_preds ) 
library(dplyr)
library(tidyr)
library(stringr)

b1_df<- test_df %>% gather( var, val , starts_with('b12')) %>% mutate( var_label = as.numeric(str_extract(pattern = '[0-9]+$', var)), m_b1_mu2 = mean(b1_mu2))

a_df<- test_df %>% gather( var, val , starts_with('a2')) %>% mutate( var_label = as.numeric(str_extract(pattern = '[0-9]+$', var)))

ggplot (b1_df, aes( x = val)) + geom_histogram() + geom_vline(aes( xintercept = m_b1_mu2)) + facet_wrap( ~ var_label, ncol = 4)

ggplot (a_df, aes( x = val)) + geom_histogram() + geom_vline(aes( xintercept = 0)) + facet_wrap( ~ var_label, ncol = 4)


testframe <- readRDS('data/temp_data/POSE_scaled_survival_dataframe.RDS')
testframe$treat_year_label
testframe %>% distinct( yid , treat_year_label) %>% arrange( yid )


hist( test_preds$mu - df$Y )
hist( test_preds$muhat - df$Yhold )
hist( test_preds$muhat2 - df$Yhold)

source( 'analysis/waic_fxns.R')

waic(testfit, llname = 'log_lik2')

summary(testfit, c('bg', 'b1_mu', 'sig_a', 'sig_b1'))$summary

test_df <- data.frame(gid = factor(df$gid), Y = df$Y, yid = df$yid, X = df$X, W.ARTR = df$W[,1], P.f.w.sp.l = df$C[, 3], `P.f.w.sp.l:logarea` = df$C[,10])
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