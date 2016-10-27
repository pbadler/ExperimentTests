rm(list = ls() )

library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 
model_table <- read.csv('output/best_WAIC_scores.csv')

model_table <- subset( model_table, vital_rate == 'recruitment' )

# log-pointwise predictive density -------------------------------------------------------# 

compute_lppd <- function( stan_fit, ll = 'log_lik' ) { 
  log_lik <- rstan::extract(stan_fit, ll)[[ll]]
  lppd <- log(colMeans(exp(log_lik)))
  lppd
} 

# ---------------------------------------------------------------------------------------------------------------------

for( do_line in 1:nrow(model_table)){ 
  
  do_model <- model_table[do_line, ]
  spp <- do_model$species
  vr <- do_model$vital_rate
  m <- do_model$model
  lambda <- do_model$lambda
  
  temp_fit <- readRDS(file = file.path( 'output/stan_fits/predictions', paste(spp, vr, m, lambda, 4, 'predict.RDS', sep = '_')))

  # log-pointwise predictive density ------------------------------------------------------------------------------------# 

  lppd1 <- compute_lppd(temp_fit)
  lppd2 <- compute_lppd(temp_fit, 'log_lik2')  
  y_out <- data.frame(species = spp, vital_rate = vr, model = m, lambda = lambda , lppd1 = lppd1 , lppd2 = lppd2 )
  
  if( do_line == 1 ) { 
    write.table( y_out, file = file.path('output', 'lppd_scores.csv'), sep = ',', row.names = FALSE, append = FALSE )
  }else { 
    write.table( y_out, file = file.path('output', 'lppd_scores.csv'), sep = ',', col.names = FALSE, row.names = FALSE, append = TRUE )
  }
}


