rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

load('analysis/figure_scripts/my_plotting_theme.Rdata')
all_files <- dir( 'output/', 'predicted_and_observed*', full.names = T)

all_preds <- lapply(all_files, read.csv)

all_preds <- do.call(rbind, all_preds)

all_preds <- all_preds %>% dplyr::select ( - X ) # clean up 

all_preds <- 
  all_preds %>% 
  mutate ( sig = (lcl>0) == (ucl>0) ) %>% 
  mutate( sig = ifelse( sig , (mu != 0), sig )) %>% 
  arrange(species, vital_rate, par, type, Treatment ) %>% 
  gather( stat, val, lcl, mu, ucl)

oe <- 
  all_preds %>% 
  filter( type == 'observed_effect') %>% 
  spread(type, sig )  %>% 
  spread( stat, val )
  
pe <- 
  all_preds %>% 
  filter( type == 'predicted_effect' ) %>% 
  spread(type, sig) %>% 
  spread(stat, val )

all_effects <- 
  merge(pe, oe, by = c('Treatment', 'par', 'species', 'vital_rate'), all.x = T, all.y = T)

all_effects <- 
  all_effects %>% 
  unite( type , observed_effect, predicted_effect, sep = '_' )  %>% 
  mutate( type = ifelse(type == 'TRUE_FALSE' | type == 'FALSE_TRUE', 'mixed', type )) %>% 
  mutate( sig_label = ifelse(type == 'TRUE_TRUE', 'significant', 'not significant'))

overall_cor <- all_effects %>% ungroup %>% summarise( cor = cor(mu.x, mu.y))
overall_cor

overall_cor <- all_effects  %>% group_by(Treatment) %>% summarise( cor = cor(mu.x, mu.y ))
overall_cor

cors <- all_effects %>% group_by(par) %>% summarise( cor = cor(mu.x, mu.y))
cors$x.pos <- c(1.5,1.5) 
cors$y.pos <- c(1.5,1.5)
cors$label <- paste0( 'r=', round( cors$cor, 2))

all_effects <- unique(all_effects)

nrow( all_effects )


p1 <- 
  ggplot(all_effects, aes( x = mu.x, y = mu.y, color = Treatment, shape = sig_label,  group = 1 )  ) + 
  geom_point(size = 2.5, alpha = 0.7) + 
  geom_smooth(method = 'lm', se = F, color = 'darkgray', linetype = 2) +
  geom_vline(aes(xintercept = 0) , linetype = 2 , alpha = 0.1) + 
  geom_hline( aes(yintercept = 0 ), linetype = 2, alpha = 0.1)  + 
  ylab( 'Observed Treatment Effect') + 
  xlab( 'Predicted Treatment Effect') + 
  scale_x_continuous(limits = c(-2, 2))  + 
  scale_y_continuous( limits = c(-2,2)) + 
  scale_color_manual(values = my_colors[3:4])  + 
  my_theme + 
  scale_shape_manual(values = c(16,8), guide = 'none') + 
  facet_grid( . ~ par  )  + 
  coord_fixed()


p1 <- p1 + geom_text( data = cors, aes( x  = x.pos, y = y.pos , label = label , color = NULL, shape = NULL), show.legend = F) 
p1


p1 <- p1 + geom_text( data = subset( all_effects, type == 'TRUE_TRUE' & species == 'POSE' & vital_rate == 'recruitment' & Treatment == 'Irrigation'), 
                aes( label = paste(species, vital_rate, sep = '\n' )) , adj = -0.1, nudge_y = -0.2, show.legend = F) 

p <- 
  ggplot( all_effects , aes( x = mu.x , y = mu.y, color = Treatment ) ) + 
  geom_point()  + 
  #geom_vline(aes(xintercept = 0) , linetype = 2 ) + 
  #geom_hline( aes(yintercept = 0 ), linetype = 2)  + 
  scale_x_continuous(limits = c(-2, 2))  + 
  scale_y_continuous( limits = c(-2,2)) + 
  xlab( 'Predicted Effect') + 
  ylab( 'Observed Effect') + 
  scale_color_manual(values = my_colors[3:4]) 

g <- all_effects %>% group_by(type ) %>% do( p = p %+% . + ggtitle(.$type) )

g$p



png( 'figures/parameter_predictions.png' , height = 6, width = 8 , res = 300 , units = 'in')

print( p1 + labs( color = '' )) 

dev.off()
