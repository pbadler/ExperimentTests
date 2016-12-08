library(stringr)
library(dplyr)
library(tidyr)
library(rstan)

rm(list = ls())

# functions --------------------------------------------------------------------------------------- # 

get_parameter_estimates <- function( fit_file, treatment_fit_file, dat_files) { 
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
  
  b2  <- rstan::extract(readRDS(fit_file), 'b2')$b2
  Chold <- dat$Chold  
  
  if ( !is.matrix( b2)) { 
    b2 <- matrix(b2)
    Chold <- matrix(dat$Chold)
    colnames(Chold) <- names( dat$Ccenter )
  }
   
  VWC_vars <- grep ('VWC', colnames(Chold))
  
  VWC <- Chold[, VWC_vars, drop = F]
  IFX_vars <- grep ('logarea.t0', colnames(VWC))
  
  b2VWC <- b2[, VWC_vars, drop = F]
  
  VWC1 <- VWC[, -IFX_vars, drop = F]
  VWC2 <- VWC[, IFX_vars, drop = F]
  b2VWC1 <- b2VWC[, -IFX_vars, drop = F]
  b2VWC2 <- b2VWC[, IFX_vars, drop = F]
  
  if (length(IFX_vars) == 0) {
    VWC1 <- VWC
    b2VWC1 <- b2VWC
  }
  
  TE1 <- matrix(NA, ncol = nrow(VWC1), nrow = nrow(b2))
  TE2 <- matrix(NA, ncol = nrow(VWC2), nrow = nrow(b2))
  
  for (j in 1:nrow(TE1)) {
    TE1[j,] <- VWC1 %*% b2VWC1[j,]
    TE2[j,] <- VWC2 %*% b2VWC2[j,]
  }
  
  if (vr != 'recruitment') {
    X <- dat$Xhold
    TE2 <-
      sweep(TE2, 2, X, '/') # divide TREATMENT EFFECT by size to isolate the climate effects
  }
  
  TE1_df <- data.frame( base_df,  parameter = 'a (intercept)', t(TE1) ) 
  TE2_df <- data.frame( base_df, parameter = 'b1 (size)', t(TE2))
  
  TE_df <- rbind(TE1_df, TE2_df)
  
  pred_df <- 
    TE_df %>% 
    mutate( type = 'predicted_effect') %>% 
    gather( iteration, val , starts_with('X'))  
  
  pred_df$Treatment <- factor(pred_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
    
  avg_pred <- 
    pred_df %>% 
    filter( year > 2010 & Group == 1) %>% 
    group_by(species, vital_rate, Treatment, parameter, type , iteration ) %>% 
    summarise( val = mean(val)) %>% 
    spread(Treatment, val ) %>% 
    mutate( Drought = Drought - Control, Irrigation  = Irrigation - Control ) %>%
    gather( Treatment, val ,  Control:Irrigation ) %>% 
    filter( Treatment !=  'Control' ) 

  # get observed treatment effects ---------------------------------------- # 
  bt <- rstan::extract(readRDS(treatment_fit_file), 'bt')$bt 
  
  if(vr != 'recruitment' ) {
    obs_effect <-  
      data.frame( species = spp, 
                  vital_rate = vr, 
                  type = 'observed_effect',
                  parameter = c('a (intercept)','a (intercept)', 'b1 (size)', 'b1 (size)') ,  
                  Treatment = c('Drought', 'Irrigation', 'Drought', 'Irrigation'), 
                  t(bt) )
  }else{ 
    obs_effect <-  
      data.frame( species = spp, 
                  vital_rate = vr, 
                  type = 'observed_effect',
                  parameter = c('a (intercept)','a (intercept)'),  
                  Treatment = c('Drought', 'Irrigation'), 
                  t(bt) )
  }
  
  obs_effect <-  
    obs_effect %>% 
    gather( iteration, val , starts_with('X'))
  
  all_effects <- bind_rows( avg_pred, obs_effect)
  
  return(all_effects)
}

# ---------------------------------------------------------------------------------------- # 

fit_files <- dir( 'output/stan_fits', '*climate_fit.RDS', full.names = T)
fit_files <- fit_files[ -grep('best', fit_files)]
dat_files <- dir( 'data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
treatment_stan_fit  <- dir( 'output/stan_fits', 'treatment_fit', full.names = T)

all_effects <- list()

for( i in 1:length(fit_files)) {   
  
  all_effects[[i]] <- get_parameter_estimates( fit_files[i], treatment_stan_fit[i], dat_files )

}

all_effects <- do.call(rbind, all_effects)

par_plot_df <- 
  all_effects %>% 
  filter( !(vital_rate == 'recruitment' & parameter == 'b1 (size)')) %>% 
  group_by( species, vital_rate, Treatment, type, parameter ) %>% 
  summarise( center = mean(val), lbci = quantile( val , 0.025), ubci = quantile(val, 0.975))

par_plot <- 
  ggplot( subset( par_plot_df, species == 'ARTR' & vital_rate == 'growth'), 
          aes( x = Treatment, y = center, ymin = lbci, ymax = ubci, color = Treatment, shape = type )) + 
  geom_point(position = position_dodge(width = 0.5), size = 3) + 
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.2) + 
  geom_hline( aes( yintercept = 0 ), linetype = 2 , alpha = 0.5) + 
  facet_grid( parameter ~ .  ) + 
  ylab ( 'Mean effect (+/- 95% Bayesian Credible Interval)') + 
  theme_bw()
  
gg <- par_plot_df %>% group_by( species, vital_rate ) %>% do( gg = par_plot %+% . + ggtitle(paste( 'Parameter estimates for', .$species, .$vital_rate, sep = ' ')))

pdf( 'figures/observed_vs_predicted_treatment_parameters.pdf', width = 8, height = 5)
  gg$gg
dev.off()


compare_plot_df <- 
  par_plot_df %>% 
  gather(stat, val, center, center:ubci) %>% 
  filter( stat == 'center')%>%
  spread(type, val)  

pdf( 'figures/compare_predicted_vs_observed_treatment_parameters.pdf', width = 8, height = 5)
  
ggplot( compare_plot_df, aes( x = predicted_effect, y = observed_effect))  + 
  geom_point() + 
  facet_wrap( ~ species ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')
  
ggplot( compare_plot_df, aes( x = predicted_effect, y = observed_effect))  + 
  geom_point() + 
  facet_wrap( ~ vital_rate ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')

ggplot( compare_plot_df, aes( x = predicted_effect, y = observed_effect))  + 
  geom_point() + 
  facet_wrap( ~ Treatment ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')

dev.off()


