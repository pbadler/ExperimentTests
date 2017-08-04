# plot obs v predicted treatment parameters together 
library(ggplot2)
library(dplyr)
library(tidyr)

fils <- dir('output/', 'predicted_and_observed_par_estimates.*.csv', full.names = T)

dat <- lapply( fils, read.csv)

dat <- do.call(rbind, dat)

load('analysis/figure_scripts/my_plotting_theme.Rdata')

dat$par <- factor(dat$par, labels = c('Intercept', 'Size x Treatment'))
dat$type <- factor(dat$type, levels = c('predicted_effect', 'observed_effect'), ordered = T)

blank_plot <- 
  ggplot( data = dat, aes( x = Treatment, y = mu, ymin = lcl, ymax = ucl, color = Treatment, shape = type)) + 
  geom_point(position = position_dodge(width = 0.3) , size = 4, color = NA) + 
  geom_errorbar(position = position_dodge(width = 0.3), width = 0.3, color = NA) + 
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.3 ) +
  scale_color_manual(values = my_colors[3:4]) + 
  facet_grid( par ~ species)  + 
  ylab('Parameter estimate (+/- 95% BCI)') +
  my_theme + theme(strip.background = element_blank())
  
png('figures/blank_treatment_effects.png', height = 5.8, width = 8 , units = 'in', res = 300)
print(blank_plot)
dev.off()

par_plot <- 
  ggplot( data = dat, aes( x = Treatment, y = mu, ymin = lcl, ymax = ucl, color = Treatment, shape = type)) + 
  geom_point(position = position_dodge(width = 0.3) , size = 4) + 
  geom_errorbar(position = position_dodge(width = 0.3), width = 0.3) + 
  geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.3 ) +
  scale_color_manual(values = my_colors[3:4]) + 
  scale_shape_manual(values = c(17,19) ) + 
  facet_grid( par ~ species)  + 
  ylab('Treatment effects (+/- 95% Bayesian Credible Interval)') +
  my_theme + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank())


p1 <- dat %>% filter( type == 'predicted_effect')  %>% group_by(vital_rate) %>% do( p = par_plot %+% . )
p1$p

p2 <- dat %>% group_by(vital_rate) %>% do( p = par_plot %+% . )

p2$p[1] 

png('figures/pred_v_obs_growth1.png', height = 5.26, width = 7.2, units = 'in', res = 300)
print(p1$p[[1]] + theme( axis.text = element_blank(), axis.title = element_blank(), strip.text = element_blank()) + scale_color_manual(values = my_colors[3:4], guide = F) +   scale_shape_manual(values = c(17,19), guide = F ) + ylim(c(-1.5, 1.75))) 
dev.off()

png('figures/pred_v_obs_recruitment1.png', height = 5.26, width = 7.2, units = 'in', res = 300)
print(p1$p[[2]] + theme( axis.text = element_blank(), axis.title = element_blank(), strip.text = element_blank()) + scale_color_manual(values = my_colors[3:4], guide = F) +   scale_shape_manual(values = c(17,19), guide = F ) )
dev.off()

png('figures/pred_v_obs_survival1.png', height = 5.26, width = 7.2, units = 'in', res = 300)
print(p1$p[[3]] + theme( axis.text = element_blank(), axis.title = element_blank(), strip.text = element_blank()) + scale_color_manual(values = my_colors[3:4], guide = F) +   scale_shape_manual(values = c(17,19), guide = F ) )
dev.off()


png('figures/pred_v_obs_growth.png', height = 5.26, width = 7.2, units = 'in', res = 300)
print(p2$p[[1]] + theme( axis.text = element_blank(), axis.title = element_blank(), strip.text = element_blank()) + scale_color_manual(values = my_colors[3:4], guide = F) +   scale_shape_manual(values = c(17,19), guide = F ) )
dev.off()

png('figures/pred_v_obs_recruitment.png', height = 3.2, width = 7.2, units = 'in', res = 300)
print(p2$p[[2]] + theme( axis.text = element_blank(), axis.title = element_blank(), strip.text = element_blank()) + scale_color_manual(values = my_colors[3:4], guide = F) +   scale_shape_manual(values = c(17,19), guide = F ) )
dev.off()

png('figures/pred_v_obs_survival.png', height = 5.26, width = 7.2, units = 'in', res = 300)
print(p2$p[[3]] + theme( axis.text = element_blank(), axis.title = element_blank(), strip.text = element_blank()) + scale_color_manual(values = my_colors[3:4], guide = F) +   scale_shape_manual(values = c(17,19), guide = F ) )
dev.off()

