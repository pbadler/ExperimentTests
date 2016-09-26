#!/usr/bin/env 

##########################################################################
#
##  Calculate WAIC for finished models on HPC  
##
# 
##########################################################################

rm(list = ls())
library(rstan)
library(matrixStats)

# function ----------------------------------------------------------------------------- # 

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

# --------------------------------------------------------------------------------------# 

args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 2){ 
  stop('####### Incorrect number of arguments supplied ####### \n
       ####### Arguments required:  
       #######  path: path for folder containing fits 
       #######  do_file:  file number, 1-f where f is the total number of RDS stan fits in the folder ')
}else if (length(args) == 2){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#
  args <- commandArgs(trailingOnly = TRUE)
  
  print(paste( c('path', 'do_file') , args))
  
  setwd(args[1])  # set to directory with the "fits"
  
  do_file <- as.numeric(eval(parse (text = args[2])))
 
}

files <- dir( '.' , pattern = '*.RDS', full.names = T)[ do_file ]

print(files)

file_df <- data.frame(fn = files )


if( nrow(file_df) == 1){
  
    file_df <-  data.frame(file_df, waic(readRDS(file = as.character( file_df$fn) )))
    saveRDS(file_df, file.path('..', 'WAIC_cache', gsub(basename( as.character(file_df$fn)), pattern = '.RDS', replacement = '_WAIC.RDS')))
    
  }else(stop('incorrect number of files provided'))


