load('data/temp_data/master_list.Rdata')
master_list$sd_vec

df <- readRDS('data/temp_data/survival_data_lists_for_stan.RDS')[['ARTR']]
colnames(df$W)
colnames(df$C)
summary(df$C[,10])

gm <- model.matrix.lm(~ factor(df$gid))

df$gm <- as.matrix( gm  ) 

library(rstan)

head(df$C)

df$tau_beta <- master_list$sd_vec[1, 2]

testfit <- stan('analysis/survival/model_survival_1.stan', chains = 1, data = df, cores = 1, iter = 2000)

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