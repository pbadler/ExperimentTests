library(ggmcmc)
library(rstan)
library(matrixStats)

waic <- function(stanfit){
  
  # from Gelman 2014
  # http://www.stat.columbia.edu/~gelman/research/unpublished/waic_stan.pdf
  
  log_lik <- extract (stanfit, "log_lik")$log_lik
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                    colMeans(loo_weights_regularized))
  
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic = waic,lppd = lpd, p_waic = p_waic,elpd_waic = elpd_waic, p_loo = p_loo, elpd_loo = elpd_loo)
  
  total <- colSums(pointwise)
  
  se <- sqrt(n*colVars(pointwise))
  names(se) <- paste(names(total), 'se', sep = '_')
  return(total = data.frame( t( total ), t(se) ))
}

do_spp <- 'ARTR'
do_vital_rate <- 'growth' 
do_model <- 2

dl <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')
inits <- readRDS('data/temp_data/growth_init_vals.RDS')

# regularization based on Gerber et al. 2015 ---------------------------------------------------------------------# 
nlambda <- 30
lambda.set <- exp(seq(-5, 15, length=nlambda))
sd_vec <- sqrt(1/lambda.set) # use sd for stan normal distribution 
# ----------------------------------------------------------------------------------------------------------------# 

dl[[do_spp]]$tau_beta <- sd_vec[1]

out1 <- rstan::stan('analysis/growth/model_growth_2.stan', data = dl[[do_spp]], init = rep(inits[[do_spp]][do_model], 4), pars = c('log_lik', 'b2') , cores = 4 )

waic_1 <- round( waic(out1), 2)
waic_1

dl[[do_spp]]$tau_beta <- sd_vec[30]

out2 <- rstan::stan('analysis/growth/model_growth_2.stan', data = dl[[do_spp]], init = rep(inits[[do_spp]][do_model], 4), pars = c('log_lik', 'b2') , cores = 4 )
waic_2 <- round( waic(out2), 2)
waic_2


b2_1 <- ggs(out1, 'b2')
b2_2 <- ggs(out2, 'b2')

b2_1$model <- paste('sd prior =', sd_prior[1], 'waic =', waic_1$waic)
b2_2$model <- paste('sd prior =', sd_prior[30], 'waic =', waic_2$waic)

both_models <- rbind(b2_1, b2_2)

climate_vars <- colnames(dl[[do_spp]]$C )

both_models$par_label <- factor( both_models$Parameter, labels = climate_vars)


pdf(file.path('figures', paste('Climate regularization plot for', do_spp, do_vital_rate, do_model, '.pdf', sep = '_')), width = 8, height = 11 )

print(
  ggplot(both_models, aes( x = value ) ) + geom_density() + facet_grid(par_label ~ model) + ggtitle(paste('Regularization effect for', do_spp, do_vital_rate, 'model' , do_model))
)

dev.off()
