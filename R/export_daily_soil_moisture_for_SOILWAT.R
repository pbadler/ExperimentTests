# clean soil moisture data 

rm(list = ls () )
library( dplyr ) 
library( tidyr )

soil <- readRDS('data/temp_data/decagon_data_corrected_values.RDS')

soil %>% filter( Treatment == 'Control') %>% mutate( year = strftime( date, '%Y')) %>% group_by(plot, year ) %>% summarise( n())

soil_export <- 
  soil %>% 
  filter( measure == 'VWC', 
          stat == 'raw', 
          !is.na(v), 
          bad_values == 0, 
          good_date == 1) %>% 
  dplyr::select(port, plot, Time, date, modified_date, measure, quad, Treatment, PrecipGroup, depth, v)

dups <- soil_export %>% 
  ungroup() %>% 
  select( plot, port, Time, date, v ) %>% 
  distinct %>% 
  group_by( plot, date, port ) %>% 
  mutate( cnt  = n()) %>% 
  filter ( cnt > 1 ) %>% 
  select( plot, port, Time, date, cnt, v ) %>% 
  arrange( plot, port , date ) %>% 
  mutate(duplicate = TRUE)

data.frame( dups ) 

soil_export %>% ungroup() %>% dplyr::select(plot, port, depth ) %>% distinct() %>% group_by(plot) %>% mutate( SMS_number = n() ) 

port_labels <- data.frame( port = unique( soil_export$port), new_label = paste0('VWC_L', 1:4))

temp <- soil_export %>% left_join(port_labels, by = 'port') %>% ungroup() 

temp <- temp %>% distinct()

temp <- temp %>% select( plot, Treatment, date, v, new_label ) %>% distinct()

temp_avg <- 
  temp %>% 
  mutate( old_date = date) %>%
  mutate( date = strftime( old_date, '%Y-%m-%d'), DOY  = strftime( old_date, '%j'), year = strftime( old_date , '%Y') ) %>% 
  group_by( plot, date, Treatment, year, DOY, new_label) %>% 
  summarise( VWC = mean(v), n = n()) %>% 
  ungroup()

temp_avg <- temp_avg %>% 
  group_by( plot, year, new_label) %>% 
  mutate( scaled = scale(VWC)) %>% 
  filter(! ( plot == '7_8_C' & VWC < -0.01 & new_label == 'VWC_L3')) # filter out bad values in plot 7_8_C

library(ggplot2)

ggplot( temp_avg %>% filter( Treatment == 'Control'), aes( x = DOY, y = VWC, color = new_label )) + 
  geom_point() + facet_grid( new_label ~ plot )

ggplot( temp_avg %>% filter( Treatment == 'Control'), aes( x = date, y = VWC, color = new_label )) + 
  geom_point() + facet_grid( plot  ~  new_label )

ggplot( temp_avg %>% ungroup() %>% filter( Treatment == 'Control' & new_label %in% c('VWC_L1' , 'VWC_L2')), aes( x = date, y = scaled, group = plot))  + 
  geom_line() + facet_wrap(~plot + new_label, ncol = 1)

ggplot ( temp_avg %>% filter( Treatment == 'Control' & plot == '7_8_C'), aes( x = date, y = VWC, group = new_label, color = new_label)) + geom_line() + 
  facet_wrap(~ new_label)

control_data <- temp_avg %>% 
  ungroup( ) %>%
  select( - scaled , - n ) %>% 
  spread(  new_label, VWC) %>% 
  filter( Treatment == 'Control') %>%
  rename( Date = date, doy = DOY) 

out_list <- split( select( control_data, Date, doy, starts_with("VWC") ), control_data$plot )

for( i in 1:length(out_list) ) { 
  fname <- paste0("USSES_", names(out_list)[i], '_SoilWater.csv')
  write.csv(out_list[[i]], file.path('data/temp_data/soil_files', fname), row.names = FALSE )
}


write.table( control_data, 'data/temp_data/soil_moisture_controls.csv', sep = ',', row.names = F)

