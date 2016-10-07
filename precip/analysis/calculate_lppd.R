rm(list = ls() )

library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 
model_table <- read.csv('output/best_WAIC_scores.csv')

model_table <- subset( model_table, vital_rate != 'recruitment' )

# log-pointwise predictive density -------------------------------------------------------# 

compute_lppd <- function( stan_fit ) { 
  log_lik <- rstan::extract(stan_fit, "log_lik")$log_lik
  lppd <- log(colMeans(exp(log_lik)))
  lppd
} 

# ---------------------------------------------------------------------------------------------------------------------

for( do_line in 1:nrow(model_table)){ 
  
  do_model <- model_table[do_line, ]
  spp <- do_model$species
  vr <- do_model$vital_rate
  m <- do_model$model
  prior <- do_model$prior
  
  temp_fit <- readRDS(file = file.path( 'output/stan_fits/predictions', paste(spp, vr, m, prior, 4, 'predict.RDS', sep = '_')))

  df1 <- readRDS(file.path( 'data/temp_data', paste0 ( spp, '_', vr, '.RDS')) )

  df <- readRDS(paste0( 'data/temp_data/', vr, '_data_lists_for_stan.RDS'))[[spp]]
  
  if( vr %in% c('growth', 'survival')){ 
    y_out <- data.frame(
      species = spp ,
      vital_rate = vr, 
      model = m , 
      prior = prior, 
      obs_id = 1:length(df$y_holdout),
      Y = df$y_holdout,
      X = df$Xhold,
      Period = 'Modern',
      treatment = factor( df$treat_out, labels = c('Control', 'Drought', 'Irrigation')),
      G = df$G,
      year = df$yid_out, 
      cyear = df$year_out, 
      trackid = df$trackid_out, 
      quad = df$quad_out)
    
    y_training <- data.frame( 
      species = spp, 
      vital_rate = vr, 
      model = m , 
      prior = prior, 
      obs_id = 1:length(df$Y), 
      Y = df$Y,
      X = df$X,
      Period = 'Historical', 
      treatment = 'Control',
      G = df$G,
      year = df$yid, 
      cyear = df$year, 
      trackid = df$trackid, 
      quad = df$quad)
  } 
    
  # log-pointwise predictive density ------------------------------------------------------------------------------------# 

  lppd <- compute_lppd(temp_fit)
  
  muhat <- rstan::extract(temp_fit, 'muhat')$muhat
  mu    <- rstan::extract(temp_fit, 'mu')$mu
  
  y_out <-  cbind ( y_out ,  t( apply (muhat, 2 ,  quantile,  c(0.025, 0.1, 0.5, 0.9, 0.975))))
  y_training <- cbind ( y_training, t(apply( mu, 2, quantile,  c(0.025, 0.1, 0.5, 0.9, 0.975))))

  y_out$lppd <- lppd  
  y_training$lppd <- NA
  
  y_out <- rbind( y_out , y_training ) 
  
  if( do_line == 1 ) { 
    write.table( y_out, file = file.path('output', 'lppd_scores.csv'), sep = ',', row.names = FALSE, append = FALSE )
  }else { 
    write.table( y_out, file = file.path('output', 'lppd_scores.csv'), sep = ',', col.names = FALSE, row.names = FALSE, append = TRUE )
  }
}


