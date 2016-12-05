library(stringr)
library(dplyr)
library(tidyr)
library(rstan)

rm(list = ls())

# functions --------------------------------------------------------------------------------------- # 

get_climate_contributions <- function( fit_file, treatment_fit_file, dat_files) { 
  bname <- basename(fit_file)
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]

  dat <- readRDS( dat_files[grep( vr, dat_files)] )[[spp]]
  
  # make dataframe for predictions ------------------------------------
  base_df <- data.frame(species = spp, 
                        vital_rate = vr,
                        year = dat$yearhold, 
                        Group = dat$gidhold, 
                        Treatment = dat$treathold, 
                        obs = 1:length(dat$treathold))
  
  CC    <- rstan::extract(readRDS(fit_file), 'climhat')$climhat

  pred_df <- 
    data.frame( base_df, t(CC) ) %>% 
    mutate( type = 'predicted_effect') %>% 
    gather( iteration, val , starts_with('X'))  
  
  # get observed treatment effects ---------------------------------------- # 

  YE <- rstan::extract(readRDS(treatment_fit_file), 'year_effect')$year_effect 
  TE <- rstan::extract(readRDS(treatment_fit_file), 'treatEff')$treatEff
  
  CC_obs <- YE + TE # add treatment and year effect to get total contribution from year and treatment 
  
  obs_df <- 
    data.frame( base_df, t(CC_obs) ) %>% 
    mutate( type = 'observed_effect') %>% 
    gather( iteration, val , starts_with('X'))  
  
  all_df <- bind_rows( pred_df , obs_df )
  
  all_df$Treatment <- factor(all_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
  
  avg_CC <- 
    all_df %>% 
    group_by(species, vital_rate, year, Treatment, type, iteration ) %>% 
    summarise( val = mean(val)) 
  
  return(avg_CC)
}

# ---------------------------------------------------------------------------------------- # 
fit_files <- dir( 'output/stan_fits', '*climate_fit.RDS', full.names = T)
dat_files <- dir( 'data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
treatment_stan_fit  <- dir( 'output/stan_fits', 'treatment_fit', full.names = T)

all_effects <- list()

for( i in 1:length(fit_files)) {   
  all_effects[[i]] <- get_climate_contributions( fit_files[i], treatment_stan_fit[i], dat_files )
}

all_effects <- do.call(rbind, all_effects)

error <- 
  all_effects %>% 
  spread( type, val ) %>% 
  mutate( error = predicted_effect - observed_effect ) %>%
  gather( type , val, observed_effect:error)

write.csv(error , 'output/climate_effect_prediction_error.csv')

avg_effects <- 
  all_effects %>% 
  group_by( species, vital_rate, year, Treatment, type) %>% 
  summarise( center = mean(val), lbci = quantile( val , 0.025), ubci = quantile(val, 0.975))

plot_CC <- 
  ggplot( subset( avg_effects, species == 'ARTR' & vital_rate == 'growth'), 
          aes( x = year, y = center, ymin = lbci, ymax = ubci, color = Treatment, shape = type )) + 
  geom_point(position = position_dodge(width = 0.5), size = 3) + 
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.2) + 
  geom_hline( aes( yintercept = 0 ), linetype = 2 , alpha = 0.5) + 
  facet_grid( Treatment  ~ . ) + 
  ylab ( 'Average climate contribution (+/- 95% Bayesian Credible Interval)') + 
  theme_bw()
  
gg <- 
  avg_effects %>% 
  group_by( species, vital_rate ) %>% 
  do( gg = plot_CC %+% . + ggtitle(paste( 'Climate effects for', .$species, .$vital_rate, sep = ' ')))

pdf( 'figures/observed_and_predicted_climate_contributions.pdf', width = 8, height = 5)
  gg$gg
dev.off()

compare_plot_df <- 
  avg_effects  %>%
  dplyr::select(species, vital_rate, year, Treatment, type, center, lbci , ubci ) %>% 
  gather( stat, val, center:ubci) %>% 
  unite( type, stat, type ) %>%
  spread(type, val )


cp <- 
  ggplot(compare_plot_df, aes( x = center_predicted_effect, 
                             y = center_observed_effect, 
                             ymin = lbci_observed_effect, 
                             ymax = ubci_observed_effect, 
                             xmin = lbci_predicted_effect, 
                             xmax = ubci_predicted_effect)) + 
  geom_point() +
  geom_errorbar(alpha = 0.1) + 
  geom_errorbarh(alpha = 0.1)


write.csv( compare_plot_df, 'output/avg_observed_and_predicted_climate_effects.csv' )


pdf( 'figures/compare_observed_vs_predicted_climate_contributions.pdf', width = 8, height = 5)
  
  cp + 
  facet_wrap( ~ species ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')
  
  cp + 
  facet_wrap( ~ vital_rate ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')

  cp + 
  facet_grid( species ~ Treatment ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')
  
  cp + 
  facet_grid( species ~ vital_rate ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')

dev.off()


