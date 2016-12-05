rm(list = ls() )

library(ggmcmc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(rstan)
library(gridExtra)

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')

dir('output/stan_fits/')
climate_models    <- dir('output/stan_fits', '.*_survival_climate_fit.RDS', full.names = TRUE)
treatment_models  <- dir('output/stan_fits', '.*_survival_treatment_fit.RDS', full.names = T)
year_models       <- dir('output/stan_fits', '.*_survival_year_effects_fit.RDS', full.names = T)

dffile <- dir('data/temp_data', 'modified_survival_data_lists_for_stan.RDS', full.names = T)


extract_predictions <- function( climate_model, treatment_model, year_model, dat ){ 
  
  bname <- basename(climate_model)
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]

  # make dataframe for predictions ------------------------------------
  base_df <- data.frame(species = spp, 
                        vital_rate = vr,
                        Y = dat$Yhold,
                        logarea.t0 = dat$Xhold,
                        year = dat$yearhold, 
                        Group = dat$gidhold, 
                        quad = dat$quadhold, 
                        trackid = dat$trackidhold, 
                        Treatment = dat$treathold, 
                        obs = 1:length(dat$treathold))
  
  # extract predictions from climate model --------------------------- 
  mu <- rstan::extract(readRDS(climate_model), 'muhat')$muhat # extract linear predictor for held out data 
  
  climate_pred <- data.frame(base_df, model = 'climate_model', t(mu))
  
  climate_pred <- 
    climate_pred %>% 
    gather( 'iteration', mu, starts_with('X'))
  
  
  # extract predictions from treatment model -------------------------
  
  mu <- rstan::extract(readRDS(treatment_model), 'mu')$mu  # extract linear predictor 
  
  treatment_pred <- data.frame(base_df, model = 'treatment_model', t(mu))
  
  treatment_pred <- 
    treatment_pred %>% 
    gather( 'iteration', mu, starts_with('X'))

  # extract predictions from year model (no treatment or climate effects )  
  
  mu <- rstan::extract(readRDS(year_model), 'muhat')$muhat  # extract linear predictor 
  
  year_pred <- data.frame(base_df, model = 'year_model', t(mu))
  
  year_pred <- 
    year_pred %>% 
    gather( 'iteration', mu, starts_with('X'))
  
  # ----------------------------------------------------------------------- # 
  
  all_pred <- do.call( rbind , list( climate_pred, treatment_pred, year_pred ))
  
  all_pred$Treatment <- factor(all_pred$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
  
  return(all_pred)
  
} 


all_df <- list()

for( i in 1:length(climate_models)){ 
  all_df[[i]] <- extract_predictions(climate_models[i] , treatment_models[i] , year_models[i], readRDS(dffile)[[i]])
}

df <- 
  do.call(rbind, all_df) %>% 
  mutate( error = mu - Y, `size error (predicted - observed)` = exp(mu) - exp(Y) )

avg_error_per_iteration <- 
  df %>% 
  group_by(Treatment, species, model, iteration ) %>% 
  summarise(mean_error = mean(error))
    
gg_it <- 
  ggplot( subset ( avg_error_per_iteration, species == 'ARTR'), aes( x = mean_error))  + 
  geom_density(alpha = 0.5) + 
  geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.5) + 
  facet_grid( model ~ . )

p1 <- avg_error_per_iteration %>% 
  group_by(species) %>% 
  do(gg = gg_it %+% .  + ggtitle (paste( 'Average error per iteration for', .$species, 'survival')))

avg_error_per_obs <- 
  df %>%
  group_by( Treatment, species, model, obs ) %>% 
  summarise( mean_error = mean(error))

gg_obs <- 
  ggplot( subset( avg_error_per_obs, species == 'ARTR'), aes ( x = mean_error) ) + 
  geom_density (alpha = 0.5) + 
  geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.5) + 
  facet_grid( model ~ . ) 

p2 <- avg_error_per_obs %>% 
  group_by(species) %>% 
  do(gg = gg_obs %+% .  + ggtitle (paste( 'Average error per observation for', .$species, 'survival')))

df$size_class <- cut(df$logarea.t0, 3, labels = c('small', 'medium', 'large') )

avg_error_by_size <- 
  df %>%
  group_by( Treatment, species, size_class, model, obs ) %>% 
  summarise( mean_error = mean(error))

gg_size <- 
  ggplot( subset( avg_error_by_size, species == 'ARTR'), aes ( x = mean_error) ) + 
  geom_density (alpha = 0.5) + 
  geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.5) + 
  facet_grid( model ~ size_class ) 

p3 <- avg_error_by_size %>% 
  group_by(species) %>% 
  do(gg = gg_size %+% .  + ggtitle (paste( 'Average error per observation for', .$species, 'survival')))

pdf(file.path( 'figures', 'predictions',  paste('error_in_survival_predictions.pdf', sep = '_')), height = 8, width = 11)  

print( 
p1$gg
)
print(
p2$gg
)
print(
p3$gg
)

dev.off()
