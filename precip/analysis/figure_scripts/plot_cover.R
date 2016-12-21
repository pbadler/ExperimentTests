# 
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

load('analysis/figure_scripts/my_plotting_theme.Rdata')

# merge observed and predicted cover ---------------------------------------------------

obs_f <- dir('output/ibm/simulations', 'observed', full.names = T)
pred_f <- dir('output/ibm/simulations', 'one_step_ahead', full.names = T)

obs <- lapply( obs_f, read.csv) 
obs <- do.call(rbind, obs)

pred <- lapply( pred_f, readRDS)
pred <- do.call(rbind, pred )

pred$year <- pred$year + 1  # predicted cover is for the following year 

obs_cover <- 
  obs %>% 
  group_by(species, year, Treatment ) %>% 
  summarise( cover_obs = mean(area) )

pred_cover <- 
  pred %>% 
  group_by(species, model, year, Treatment, simulation) %>%
  summarise( cover_avg = mean(area)) %>% 
  group_by( species, model, year, Treatment) %>%
  summarise( ucl = quantile( cover_avg, 0.75), lcl = quantile( cover_avg, 0.25), median = median(cover_avg), mean = mean(cover_avg))

added <- expand.grid( year = 2011 , species = unique(pred_cover$species), model = unique(pred_cover$model), Treatment = unique(pred_cover$Treatment))

pred_cover <- merge( added, pred_cover , all.x = T, all.y = T)

all_cover <- merge( obs_cover, pred_cover, all.x = T)

all_cover <- 
  all_cover %>% 
  gather( predictions, val, ucl:mean ) %>% 
  mutate(val = ifelse( year == 2011, cover_obs, val ) ) %>% 
  spread(predictions, val)

all_cover <- 
  all_cover %>% gather( type, val, cover_obs, mean, median )

all_cover <- all_cover %>% filter( year > 2010 )

all_cover <- subset( all_cover , type != 'median' )

all_cover$type  <- factor( all_cover$type, labels = c('observed', 'predicted'))

all_cover$val <- 100*(all_cover$val/10000)

my_plot <- 
  ggplot( subset( all_cover, model == 'climate'), aes( x = year, y = val, linetype = type, color = Treatment)) + 
  geom_point() + 
  geom_line() +
  ylab( 'Average cover (%)') + 
  scale_color_manual(values = my_colors[2:4])  + 
  labs(linetype = NULL) +
  my_theme + theme(legend.position = c(0.8,0.8) , legend.background = element_rect(colour = "gray"))


g <- all_cover %>% filter( model == 'climate') %>% group_by(species) %>% do( p = my_plot %+% . + ggtitle (paste0(.$species) ) )

g$p

png( 'figures/predicted_and_observed_cover.png', height = 7, width = 12, units = 'in', res = 300)
print( 
grid.arrange( g$p[[1]] + guides( color = 'none', linetype = 'none') + ylim(0, 21), 
              g$p[[2]] + guides( color = 'none', linetype = 'none') + scale_y_continuous(name = NULL, limits = c(0,3.2)),  
              g$p[[3]] + guides( color = 'none', linetype = 'none') + ylim( 0, 3.2) + scale_y_continuous(name = NULL, limits = c(0,3.2)), 
              g$p[[4]] + scale_y_continuous(name = NULL, limits = c(0,3.2)), 
              nrow = 1)
)
dev.off()


print( 
  grid.arrange( g$p[[1]] + guides( color = 'none', linetype = 'none') + ylim(0, 21) + scale_x_continuous(name = NULL), 
                g$p[[2]] + guides( color = 'none', linetype = 'none') + ylim( 0, 3.2) + scale_x_continuous(name = NULL),  
                g$p[[3]] + guides( color = 'none', linetype = 'none') + ylim( 0, 3.2) + scale_x_continuous(name = NULL), 
                g$p[[4]] + ylim( 0, 3.2), 
                ncol = 1)
)
