
waic <- function(stanfit, llname = 'log_lik'){
  
  require(rstan)
  require(matrixStats)
  
  # from Gelman 2014
  # http://www.stat.columbia.edu/~gelman/research/unpublished/waic_stan.pdf
  
  log_lik <- rstan::extract (stanfit, llname)[[llname]]
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
