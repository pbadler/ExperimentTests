
###########################################################################
#
# Output monthly summaries of climate for historical Period, long-term average 
# and contemporary Period. 
#
###########################################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(gridExtra)

monthly <- readRDS('data/temp_data/monthly_climate.RDS')
annual <- readRDS('data/temp_data/annual_climate.RDS')

annual <- 
  annual %>% 
  ungroup( ) %>% 
  filter( year > 1925, year < 2016 ) %>% 
  arrange( year ) %>% 
  select( - n ) 

monthly <- 
  monthly %>% 
  ungroup() %>% 
  filter( year > 1925 & year < 2016 ) %>%
  arrange( year, month )

# plot parameters -------------------------------------------------------------------------

my_colors <- c( '#8da0cb', 'gray','#fc8d62')


# Summaries ------------------------------------------------------------------------- 

monthly_stats <- 
  monthly %>% 
  gather( var, val , TPCP:MMNT ) %>% 
  group_by( var, Period, month  ) %>% 
  summarise( avg = mean(val, na.rm = TRUE), 
             sd = sd(val, na.rm = TRUE), 
             n = n() , 
             UCL = avg + 1.96*(sd/(sqrt(n))), 
             LCL = avg - 1.96*(sd/sqrt(n)), 
             uq = quantile(val, 0.975, na.rm = TRUE ), 
             lq = quantile(val, 0.025, na.rm = TRUE),
             start_year = min(year), 
             end_year = max(year), 
             min = min(val), 
             max = max(val),
             extreme_high_year = year[ which.max( val ) ], 
             extreme_low_year = year[ which.min(val) ]) %>% 
  mutate( month_label = factor(month, labels = c(month.abb)), 
          Period_label = factor( Period, levels = c('Historical', 'not monitored', 'Modern'), ordered = TRUE), 
          var_label = factor( var , 
                              levels = c('MMNT', 'MNTM', 'MMXT', 'TPCP'), 
                              labels = c('Mean~monthly~min~temperature~degree*C', 
                                         'Mean~monthly~temperature~degree*C',
                                         'Mean~monthly~max~temperature~degree*C', 
                                         'Total~monthly~precip.~(mm)' )))

annual_stats  <- 
  annual %>% 
  mutate( TPPT = TPPT*10) %>% 
  gather( var, val, TPPT:MAT ) %>% 
  group_by( var, Period ) %>% 
  summarise( avg = mean(val), 
             sd = sd(val), 
             n = n(), 
             UCL = avg + 1.96*(sd/(sqrt(n))), 
             LCL = avg - 1.96*(sd/sqrt(n)), 
             uq = quantile( val, 0.95), 
             lq = quantile( val, 0.05), 
             start_year = min(year), 
             end_year = max(year), 
             min  = min (val ) , 
             max = max(val ), 
             extreme_low_year = year[ which.min(val ) ], 
             extreme_high_year = year[ which.max( val ) ]) %>% 
  mutate(Period = factor( Period, levels = c('Historical', 'not monitored', 'Modern'), labels = c('Historical', 'not monitored', 'Modern'), ordered = TRUE), 
         Period_lab = factor( Period, labels = c('Hist.', 'Not Mon.', 'Mod.')), 
         var_label = factor( var , levels = c('MAT', 'TPPT'), labels = c('Mean~annual~temperature~degree*C', 'Total~annual~precipitation~(mm)')))


# --------- make plots ----------------------------------------------------------------------------------------------------------------

monthly_stats_temp <- 
  monthly_stats  %>% 
  gather ( stat , val , avg:extreme_low_year) %>% 
  unite(stat, Period_label, stat ) %>%
  ungroup() %>% 
  select ( -Period) %>% 
  spread(stat, val )

monthly_long <- 
  monthly %>% 
  arrange( year, month )  %>% 
  mutate( date = as.Date(paste(year, month, '01', sep = '-'), format = '%Y-%m-%d')) %>%
  gather( var, observed, TPCP:MMNT) %>%
  left_join(monthly_stats_temp, by = c('var', 'month')) 

annual_stats_temp <- 
  annual_stats %>% 
  gather( stat, val, avg:extreme_high_year) %>% 
  unite(stat, Period, stat) %>% 
  select(- Period_lab) %>% 
  spread(stat, val)


annual_long <-
  annual %>% 
  mutate( TPPT = TPPT*10 ) %>% 
  gather( var, observed, TPPT:MAT) %>%
  left_join(annual_stats_temp, by = c('var')) 


month_ts_plot <- function( df ) { 
  
  df <- df %>% gather( stat, val , observed, Historical_avg )
  ggplot( df , aes( x = date, y = val, color = stat )) + 
    geom_point(alpha = 0.5) +
    geom_line( alpha = 0.5) + 
    scale_color_manual(values = my_colors[c(1,3)] ) + 
    labs( x = 'Date', 
          y = parse(text = as.character(df$var_label[1])))

}

monthly_ts_ppt <- month_ts_plot ( monthly_long %>% filter( Period == 'Modern' , var == 'TPCP'))
monthly_ts_meanT <- month_ts_plot (monthly_long %>% filter( Period == 'Modern', var == 'MNTM'))

month_ts_plot (monthly_long %>% filter( Period == 'Modern', var == 'MMNT'))
month_ts_plot (monthly_long %>% filter( Period == 'Modern', var == 'MMXT'))

annual_ts_plot <- function( df ) { 
  
  #df <- df %>% gather( stat, val , raw, Historical_avg, Historical_lq, Historical_uq )
  ggplot( df , aes( x = year, y = observed )) + 
    geom_hline( aes(yintercept = Historical_avg), col = my_colors[1], size = 1.2) + 
    geom_hline( aes(yintercept = Historical_lq), col = my_colors[1], lty = 2, size = 1.2) + 
    geom_hline( aes(yintercept = Historical_uq), col = my_colors[1], lty = 2, size = 1.2) +
    geom_line( ) + 
    geom_point( ) + 
    scale_color_manual(values = my_colors ) + 
    scale_y_continuous(limits = c(0, 1.2*max(df$Historical_uq, df$observed))) + 
    labs( x = 'Year', 
          y = parse(text = as.character(df$var_label[1])))
}

annual_ts_ppt <- annual_ts_plot(annual_long %>% filter( Period == 'Modern', var == 'TPPT'))

annual_ts_MAT <- annual_ts_plot(annual_long %>% filter( Period == 'Modern', var == 'MAT'))


# plot avgs ------------------------------------------------------------------

month_plot <- 
  function(df ) { 
    ggplot( df, aes( x = month_label, group = Period_label, color = Period_label,  fill = Period_label, y = avg, ymin = LCL, ymax = UCL  ) ) + 
      geom_point ( position = position_dodge(width = 0.8), size = 2) +
      geom_errorbar( position = position_dodge(width = 0.8), size = 1 ) + 
      scale_color_manual(values = my_colors) + 
      labs( x = 'Month', 
            y = parse(text = as.character( df$var_label[1] )) ) 
  }


ann_plot <-function( df ){  
  ggplot(df, aes( x = Period_lab, color = Period, y = avg, ymin = LCL, ymax = UCL )) + 
  geom_point(size = 4) + 
  geom_errorbar(size = 1.5) + 
  geom_point(aes( x = Period_lab, y = max), size = 3) + 
  geom_point(aes( x = Period_lab, y = min), size = 3) +
  geom_text( aes( x = Period_lab, y = min, label = extreme_low_year),  size = 3, vjust = -1.2, hjust = -0.1) + 
  geom_text( aes( x = Period_lab, y = max, label = extreme_high_year), size = 3, vjust = -1.2, hjust = -0.1) + 
  scale_color_manual(values = my_colors) + 
  scale_y_continuous( limits = c(min(df$min), 1.05*max(df$max))) + 
  labs( x = 'Period', 
        y = parse(text = as.character( df$var_label[1] )) ) 
}


p1 <- 
  month_plot ( subset(monthly_stats, var == 'MMNT') )  + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14, margin = margin(0, 20, 0, 0)), 
        plot.margin = margin(10, 10, 5, 20)) + 
  guides(color = 'none', fill = 'none')


p2 <- 
  ann_plot(subset( annual_stats, var == 'MAT')) + 
  theme(legend.key.width = unit(50, units = 'pt'),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, margin = margin(0, 20, 0, 0)), 
        plot.margin = margin(10,20,5, 30)) 

p3 <- 
  month_plot(subset( monthly_stats, var == 'TPCP')) + 
  guides(color = 'none', fill = 'none') + 
  theme(axis.text = element_text(size = 10), 
        axis.title.x = element_text(size = 14, margin = margin(20, 0, 0, 0)), 
        axis.title.y = element_text(size = 14, margin = margin(0, 20, 0, 0)), 
        plot.margin = margin(5, 10, 20, 20))


p4 <- 
  ann_plot( subset( annual_stats, var == 'TPPT') ) + 
  guides(color ='none', fill  = 'none') + 
  theme(axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 14, margin = margin(20, 0, 0, 0)),
        axis.title.y = element_text(size = 14, margin = margin(0, 20, 0, 0)), 
  plot.margin = margin(5, 155, 20, 25))


plot_all <- grid.arrange(arrangeGrob(p1, p2, p3, p4, ncol = 2, widths = c(0.55, 0.45), heights = c(0.45, 0.55)))

# -----output--------------------------------------------------------------------------- # 

annual_out <- 
  annual_stats %>%
  ungroup(.)  %>% 
  select(var_label, Period, avg, sd, n, UCL, LCL, start_year, end_year, min, max, extreme_low_year, extreme_high_year)

monthly_out <- 
  monthly_stats %>% 
  ungroup(.) %>% 
  select( var_label, Period, month, month_label, avg, sd, n, UCL, LCL, start_year, end_year, min, max, extreme_low_year, extreme_high_year )


# save -----------------------------------------------------------------------------------

write.csv(x =  annual_out , file = 'output/annual_climate_stats.csv' ) 

write.csv(x = monthly_out, file = 'output/monthly_climate_stats.csv')

write.csv(x = monthly_long, file = 'output/monthly_climate_with_avg_stats.csv')

write.csv(x = annual_long, file =  'output/annual_climate_with_avg_stats.csv')


ggsave(filename = 'figures/plot_climate_averages_by_month.png',device = 'png', plot = plot_all, width = 11, height = 8, dpi = 300)

ggsave(filename = 'figures/plot_annual_mean_temp_ts.png', 
       device = 'png', 
       plot = annual_ts_MAT + 
         theme( axis.text = element_text(size = 16), axis.title = element_text(size = 20)), 
       width = 8, height = 6, dpi = 300)

ggsave( filename = 'figures/plot_annual_ppt_ts.png', 
        device = 'png', 
        plot = annual_ts_ppt, width = 8, height  = 6, dpi = 300 )

ggsave(filename = 'figures/plot_monthly_mean_temp_ts.png', device = 'png', plot = monthly_ts_meanT, width = 8, height = 6, dpi = 300)

ggsave(filename = 'figures/plot_monthly_ppt_ts.png', device = 'png', plot = monthly_ts_ppt, width = 8, height = 6, dpi = 300)


# ---------------------------------------------------------------------------------------
