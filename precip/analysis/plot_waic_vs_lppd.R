#
# plot predictive accuracy with hold out data against waic of original data 
# 

rm(list = ls())
library(ggplot2)
library(dplyr)

df <- readRDS('output/prediction_lppds.RDS')

aggreement_plot <- ggplot( df , aes( x = waic , y = prediction_lppd )) + 
  geom_point() + 
  geom_text(aes(label = model , y = prediction_lppd), hjust = 2)

plots <- df %>% 
  group_by(species , vital_rate ) %>% 
  do( gg = aggreement_plot %+% . + ggtitle( paste( .$species, .$vital_rate) ))



pdf('figures/waic_score_vs_lppd.pdf', height = 8, width = 8 )

print( plots$gg ) 

dev.off()
