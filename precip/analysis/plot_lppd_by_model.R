#
# plot predictive accuracy with hold out data against waic of original data 
# 

rm(list = ls())
library(ggplot2)
library(dplyr)

model_scores <- read.csv('output/model_scores.csv')

aggreement_plot <- ggplot( model_scores, aes( x = waic , y = lppd_out )) + 
  geom_point() + 
  geom_text(aes(label = model , y = lppd_out), hjust = 2)

plots <- 
  model_scores %>% 
  group_by(species , vital_rate ) %>% 
  do( gg = aggreement_plot %+% . + ggtitle( paste( .$species, .$vital_rate) ))

pdf('figures/waic_score_vs_lppd.pdf', height = 8, width = 8 )

print( plots$gg ) 

dev.off()

all_models <- 
  ggplot(model_scores, aes( x = factor( model), y = rank_lppd)) + 
  geom_point() + 
  ylab( 'relative predictive score')

png( 'figures/plot_relative_model_performance.png', width = 8, height = 8, units = 'in', res = 300)
print( 
  all_models
)
dev.off()

png( 'figures/plot_relative_model_performance_by_species.png', width = 8, height = 8, units = 'in', res = 300)
print( 
  all_models + facet_wrap(~ species , ncol = 1 )
  )
dev.off()

png( 'figures/plot_relative_model_performance_by_vital_rate.png', width = 8, height = 8, units = 'in', res = 300)
print( 
  all_models + facet_wrap(~ vital_rate, ncol = 1 ) 
  )
dev.off()


