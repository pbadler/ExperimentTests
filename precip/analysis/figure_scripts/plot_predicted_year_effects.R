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
   
  IFX_vars <- grep ('logarea.t0', colnames(Chold))
  
  C1  <- Chold[, -IFX_vars, drop = F]
  C2  <- Chold[, IFX_vars, drop = F]
  a   <- b2[, -IFX_vars, drop = F]
  b1  <- b2[, IFX_vars, drop = F]
  
  if (length(IFX_vars) == 0) {
    C1 <- Chold
    a  <- b2
  }
  
  YE1 <- matrix(NA, ncol = nrow(C1), nrow = nrow(b2))
  YE2 <- matrix(NA, ncol = nrow(C2), nrow = nrow(b2))
  
  for (j in 1:nrow(YE1)) {
    YE1[j,] <- C1 %*% a[j,]
    YE2[j,] <- C2 %*% b1[j,]
  }
  
  if (vr != 'recruitment') {
    X <- dat$Xhold
    YE2 <-
      sweep(YE2, 2, X, '/') # divide individual TREATMENT X SIZE EFFECTS by size to isolate the slope parameter 
  }
  
  YE1_df <- data.frame( base_df,  parameter = 'a (intercept)', t(YE1) ) 
  YE2_df <- data.frame( base_df,  parameter = 'b1 (size)', t(YE2))
  
  YE_df <- rbind(YE1_df, YE2_df)
  
  pred_df <- 
    YE_df %>% 
    mutate( type = 'predicted_effect') %>% 
    gather( iteration, val , starts_with('X'))  
  
  pred_df$Treatment <- factor(pred_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
    
  avg_pred <- 
    pred_df %>% 
    group_by(species, vital_rate, year, Treatment, parameter, type, iteration ) %>% 
    summarise( val = mean(val)) 

  # get observed treatment effects ---------------------------------------- # 

  bt <- rstan::extract(readRDS(treatment_fit_file), 'bt')$bt 
  a  <- rstan::extract(readRDS(treatment_fit_file), 'a')$a
  
  df1 <- data.frame( unique( base_df[ , c('species', 'vital_rate', 'year')]), parameter = 'a (intercept)',  Treatment = 'Control' , t(a) )

  if(vr != 'recruitment' ) {
    b1 <- rstan::extract(readRDS(treatment_fit_file), 'b1')$b1
    df2 <- data.frame( unique( base_df[ , c('species', 'vital_rate', 'year')]), parameter = 'b1 (size)', Treatment = 'Control', t(b1))
    
    year_df <- rbind(df1, df2)
    
    year_df <- 
      year_df %>% 
      gather(iteration, val, starts_with('X')) %>%
      spread(Treatment, val)
     
    treatment_pars <-  
      data.frame( species = spp, 
                  vital_rate = vr, 
                  parameter = c('a (intercept)','a (intercept)', 'b1 (size)', 'b1 (size)') ,  
                  Treatment = c('Drought', 'Irrigation', 'Drought', 'Irrigation'), 
                  t(bt))
    

  }else{ 
    year_df <- df1 
    
    year_df <- 
      year_df %>% 
      gather(iteration, val, starts_with('X')) %>%
      spread(Treatment, val)
    
    treatment_pars <-  
      data.frame( species = spp, 
                  vital_rate = vr, 
                  parameter = c('a (intercept)','a (intercept)'),  
                  Treatment = c('Drought', 'Irrigation'), 
                  t(bt) )
  }
  
  treatment_pars <- 
    treatment_pars %>% 
    gather( iteration, val, starts_with('X')) %>%
    spread( Treatment, val )
  
  obs_effect <- merge( year_df , treatment_pars, by = c('species', 'vital_rate', 'parameter', 'iteration'))
  
  obs_effect <- 
    obs_effect %>% 
    mutate( type = 'observed_effect') %>%
    mutate( Drought = Drought + Control , Irrigation = Irrigation + Control ) %>%
    gather( Treatment, val, Control:Irrigation)
  
  all_effects <- bind_rows( avg_pred, obs_effect)
  
  return(all_effects)
}

# ---------------------------------------------------------------------------------------- # 

fit_files <- dir( 'output/stan_fits', '*climate_fit.RDS', full.names = T)
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
  filter( !(Treatment != 'Control' & year < 2011 )) %>% 
  group_by( species, vital_rate, year, Treatment, type, parameter ) %>% 
  summarise( center = mean(val), lbci = quantile( val , 0.025), ubci = quantile(val, 0.975))

par_plot <- 
  ggplot( subset( par_plot_df, species == 'ARTR' & vital_rate == 'growth'), 
          aes( x = Treatment, y = center, ymin = lbci, ymax = ubci, color = Treatment, shape = type )) + 
  geom_point(position = position_dodge(width = 0.5), size = 3) + 
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.2) + 
  geom_hline( aes( yintercept = 0 ), linetype = 2 , alpha = 0.5) + 
  facet_grid( year  ~ parameter   ) + 
  ylab ( 'Mean effect (+/- 95% Bayesian Credible Interval)') + 
  theme_bw()
  
gg <- par_plot_df %>% group_by( species, vital_rate ) %>% do( gg = par_plot %+% . + ggtitle(paste( 'Parameter estimates for', .$species, .$vital_rate, sep = ' ')))

pdf( 'figures/observed_vs_predicted_year_parameters.pdf', width = 8, height = 5)
  gg$gg
dev.off()

compare_plot_df <- 
  par_plot_df %>% 
  gather(stat, val, center, center:ubci) %>% 
  filter( stat == 'center')%>%
  spread(type, val)  

pdf( 'figures/compare_predicted_vs_observed_year_parameters.pdf', width = 8, height = 5)
  
ggplot( compare_plot_df, aes( x = predicted_effect, y = observed_effect))  + 
  geom_point() + 
  facet_wrap(  ~ species ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')
  
ggplot( compare_plot_df, aes( x = predicted_effect, y = observed_effect))  + 
  geom_point() + 
  facet_grid( parameter ~  vital_rate ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')

ggplot( compare_plot_df, aes( x = predicted_effect, y = observed_effect))  + 
  geom_point() + 
  facet_grid( parameter ~ Treatment ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')

ggplot( compare_plot_df, aes( x = predicted_effect, y = observed_effect))  + 
  geom_point() + 
  facet_grid( Treatment ~ year ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')

ggplot( compare_plot_df, aes( x = predicted_effect, y = observed_effect))  + 
  geom_point() + 
  facet_grid( parameter ~ year ) + 
  geom_smooth(aes(group = 1), se = F, method = 'lm')

dev.off()

