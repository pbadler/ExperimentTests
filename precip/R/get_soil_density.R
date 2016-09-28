rm(list = ls())

soil_data <- read.csv('data/exclosure_soil_samples.csv')

soil_data$`density g per cm3` <- (soil_data$dry_weight - soil_data$bag_weight)/soil_data$soil_volume_cm3

out <- data.frame( aggregate(data = soil_data, `density g per cm3` ~ depth, 'mean'  ) ) 


out$depth <- c('shallow', 'deep')
out$depth_description <- c('0-15 cm', '15-30 cm')

write.csv( out, 'data/soil_density.csv', row.names = FALSE)


