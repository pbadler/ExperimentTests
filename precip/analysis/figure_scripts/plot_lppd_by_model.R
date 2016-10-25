#
# plot predictive accuracy with hold out data against waic of original data 
# 

rm(list = ls())
library(ggplot2)
library(dplyr)

size <- 4

model_scores <- read.csv('output/model_scores.csv')


aggreement_plot <- ggplot( model_scores, aes( x = waic , y = lppd_out )) + 
  geom_point() + 
  geom_text(aes(label = model , y = lppd_out), hjust = 2)

plots <- 
  model_scores %>% 
  group_by(species , vital_rate ) %>% 
  do( gg = aggreement_plot %+% . + ggtitle( paste( .$species, .$vital_rate) ))


pdf('figures/waic_score_vs_lppd.pdf', height = size, width = size )

print( plots$gg ) 

dev.off()

all_models <- 
  ggplot(model_scores, aes( x = model, y = rank_lppd, color = species)) + 
  geom_point() + 
  ylab( 'relative predictive score')

png( 'figures/plot_relative_model_performance.png', width = 4, height = 4, units = 'in', res = 300)
print( 
  all_models
)
dev.off()

png( 'figures/plot_relative_model_performance_by_species.png', width = size, height = size, units = 'in', res = 300)
print( 
  all_models + facet_wrap(~ species , ncol = 1 )
  )
dev.off()

png( 'figures/plot_relative_model_performance_by_vital_rate.png', width = size, height = size, units = 'in', res = 300)
print( 
  all_models + facet_wrap(~ vital_rate, ncol = 1 ) 
  )
dev.off()

# get num parameters ------------------------------------------------------------------------------------------------

example <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')[['ARTR']]

model_scores <- merge( model_scores, data.frame( model = c(1:5), k =  c(0, ncol( example$C), 1, 1 + ncol(example$C), ncol(example$C) + ncol(example$W) )), by = 'model')

model_score_by_parameters <- 
  ggplot( model_scores, aes( x = k, y = rank_lppd , color = species )) + 
  geom_point() + 
  ylab( 'relative predictive score') + 
  xlab( 'Fixed climate and competition parameters')

png( 'figures/plot_relative_model_performance_by_paremeter_number.png', width = size, height = size, units = 'in', res = 300 ) 
print( 
  model_score_by_parameters 
  )
dev.off()


