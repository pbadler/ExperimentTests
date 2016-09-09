
###########################################################################
#
# Output monthly summaries of climate for historical period, long-term average 
# and contemporary period. 
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
  group_by( var, period, month  ) %>% 
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
          period_label = factor( period, levels = c('historical', 'not monitored', 'contemporary'), ordered = TRUE), 
          var_label = factor( var , 
                              levels = c('MMNT', 'MNTM', 'MMXT', 'TPCP'), 
                              labels = c('Mean~monthly~temperature~degree*C', 
                                         'Mean~monthly~max~temperature~degree*C',
                                         'Mean~monthly~min~temperature~degree*C', 
                                         'Total~monthly~precip.~(mm)' )))

annual_stats  <- 
  annual %>% 
  mutate( TPPT = TPPT*10) %>% 
  gather( var, val, TPPT:MAT ) %>% 
  group_by( var, period ) %>% 
  summarise( avg = mean(val), 
             sd = sd(val), 
             n = n(), 
             UCL = avg + 1.96*(sd/(sqrt(n))), 
             LCL = avg - 1.96*(sd/sqrt(n)), 
             uq = quantile( val, 0.975), 
             lq = quantile( val, 0.025), 
             start_year = min(year), 
             end_year = max(year), 
             min  = min (val ) , 
             max = max(val ), 
             extreme_low_year = year[ which.min(val ) ], 
             extreme_high_year = year[ which.max( val ) ]) %>% 
  mutate(Period = factor( period, levels = c('historical', 'not monitored', 'contemporary'), labels = c('Historical', 'Not monitored', 'Contemporary'), ordered = TRUE), 
         Period_lab = factor( Period, labels = c('Hist', 'Not mon.', 'Cont')), 
         var_label = factor( var , levels = c('MAT', 'TPPT'), labels = c('Mean~annual~temperature~degree*C', 'Mean~annual~precip~(mm)')))


# --------- make plots ----------------------------------------------------------------------------------------------------------------

month_plot <- 
  function(df ) { 
    ggplot( df, aes( x = month_label, group = period_label, color = period_label,  fill = period_label, y = avg, ymin = LCL, ymax = UCL  ) ) + 
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
  geom_text( aes( x = Period_lab, y = min, label = extreme_low_year), vjust = 0, hjust = -0.3) + 
  geom_text( aes( x = Period_lab, y = max, label = extreme_high_year), vjust = 0, hjust = -0.3) + 
  scale_color_manual(values = my_colors) + 
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
  select(var_label, period, avg, sd, n, UCL, LCL, start_year, end_year, min, max, extreme_low_year, extreme_high_year)

monthly_out <- 
  monthly_stats %>% 
  ungroup(.) %>% 
  select( var_label, period, month, month_label, avg, sd, n, UCL, LCL, start_year, end_year, min, max, extreme_low_year, extreme_high_year )


# save --------- 

write.csv(x =  annual_out , file = 'output/annual_climate_stats.csv' ) 

write.csv(x = monthly_out, file = 'output/monthly_climate_stats.csv')

ggsave(filename = 'figures/plot_climate_averages_by_month.png',device = 'png', plot = plot_all, width = 11, height = 8, dpi = 300)




# ---------------------------------------------------------------------------------------

MAT <- ts( ( climate %>% filter( treatment == 'control') %>% arrange( year) %>% distinct(year) %>% filter(year < 2016))$MAT )

acf(MAT)
pacf( MAT)

arima(MAT, order = c(0, 0, 0))
arima(MAT, order = c(1, 0, 0))
arima(MAT, order = c(1, 0, 1))
arima(MAT, order = c(2, 0, 0))
arima(MAT, order = c(3, 0, 0))

arima(MAT, order = c(0, 0, 1))
arima(MAT, order = c(0, 0, 2))
arima(MAT, order = c(1, 0, 3))
?arima
m <- arima(MAT, order = c(1,0,1))
plot( m$residuals)
plot( MAT, type = 'l')

t <- 1:length(MAT)
summary(lm(MAT ~ t ))
summary(lm( m$residuals ~ t ))

plot ( MAT, type = 'l', x = 1924 + t)

# 

MAP <- ts( ( climate %>% arrange( year) %>% distinct(year) %>% filter(year < 2016))$TPPT )

acf(MAP)
pacf( MAP)

arima(MAP, order = c(1, 0, 1))
arima(MAP, order = c(1, 0, 0))
arima(MAP, order = c(0, 0, 1))
arima(MAP, order = c(0, 0, 0))
m<- arima(MAP, order = c(1, 0, 0))

plot(m$residuals)
t <- 1:length(MAP)

summary(lm(residuals(m) ~ t ))
summary(lm(MAP ~ t))
