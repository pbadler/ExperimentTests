rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

files <- dir( 'output', 'climate_model_parameters.*.csv', full.names = T)
load('analysis/figure_scripts/my_plotting_theme.Rdata')

dat <- lapply( files, read.csv)

dat <- do.call(rbind, dat )

clim <- dat[ str_detect(dat$X, 'b2'),  ] 

clim$significant <- ifelse(  ( clim$X2.5. > 0 & clim$X97.5. > 0 ) | (clim$X2.5. < 0 & clim$X97.5. < 0 ), '*', ' ' )

par_names <- read.csv('output/selected_climate_covariates.csv')

df <- list()
for( i in 1:nrow(par_names)){ 
  species = par_names$species[i]
  vital_rate = par_names$vital_rate[i]
  temp_pars <- str_split( par_names$covars[i], ',' )[[1]]
  X = paste0('b2[', 1:length(temp_pars), ']')
  df[[i]] <- data.frame(X = X, par_name = temp_pars, species = species, vital_rate = vital_rate )
}

all_par_names <- do.call(rbind, df)
clim <- merge( clim, all_par_names)

clim <- clim %>% arrange( vital_rate, species, X)

clim$lbcl95 <- clim$X2.5. 
clim$ubcl95 <- clim$X97.5.

clim$type <- NA
clim$type[ grep( 'x', clim$par_name ) ] <- 'climate x size'
clim$type[ -grep( 'x', clim$par_name ) ] <- 'climate intercept'

clim <- clim[, c('species', 'vital_rate', 'type', 'par_name', 'mean', 'se_mean', 'lbcl95', 'ubcl95', 'significant')]


clim$climate <- str_extract(clim$par_name, 'VWC\\.[a-z]+\\.[l0-1]+')

varlist <- expand.grid( pre = 'VWC', season = c('sp', 'su', 'f', 'w'), lag = c('l', 0, 1), species = c('ARTR', 'HECO', 'POSE', 'PSSP'), vital_rate = c('growth', 'recruitment', 'survival') )

varlist <- varlist %>% unite(climate, pre, season, lag, sep = '.',remove = F )

varlist <- varlist %>% arrange( lag )

clim <- merge( varlist, clim , all.x = T)

clim$climate <- str_replace( clim$climate, 'VWC.', '') 
clim$climate <- str_replace( clim$climate, 'l', 'lag')

clim$climate <- factor( clim$climate, levels = c('su.lag', 'f.lag', 'w.lag', 'sp.lag', 
                                 'su.0', 'f.0', 'w.0', 'sp.0', 
                                 'su.1', 'f.1', 'w.1', 'sp.1' ), ordered = T )


clim <- subset(clim , season != 'w')
#clim <- clim[complete.cases(clim), ]

effect_plot <- 
  ggplot( subset( clim, vital_rate == 'growth'), aes( x = climate, y = mean , ymin = lbcl95, ymax = ubcl95, color = type)) + 
  geom_point(position = position_dodge(width =0.5)) + 
  geom_errorbar(position = position_dodge(width =0.5), width = 0.5) + 
  geom_hline( aes( yintercept = 0 ), linetype = 2 , alpha = 0.5) + 
  facet_grid( species ~ .  ) + 
  ylab ( 'Mean effect (+/- 95% Bayesian Credible Interval)') + 
  xlab( '') + 
  labs(color = '') + 
  scale_color_manual(values = c(1,'red')) + 
  my_theme

gg <- clim %>% group_by(vital_rate) %>% do(gg =  effect_plot %+% .  + ggtitle( paste0('Effects on ', .$vital_rate)))

png( 'figures/climate_effect_growth.png', height = 5, width = 7, res = 300, units = 'in' ) 
print( gg$gg[[1]] ) 
dev.off() 
png( 'figures/climate_effect_recruitment.png', height = 5, width = 7, res = 300, units = 'in' ) 
print( gg$gg[[2]] + guides(color = 'none'))  
dev.off() 
png( 'figures/climate_effect_survival.png', height = 5, width = 7, res = 300, units = 'in' ) 
print( gg$gg[[3]] ) 
dev.off() 

