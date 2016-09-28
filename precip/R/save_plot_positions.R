pts <- read.csv('data/USSES_garmin_points.csv')


qinfo <- read.csv('~/driversdata/data/idaho_modern/quad_info.csv')


head(qinfo)
head(pts)

pts$name

Xplots <- subset(qinfo, Treatment != 'Control')

Xpts <- pts[ grep ( pts$name , pattern = '^X'), ] 

head(Xpts)

Xpts <- Xpts[ , c('lat', 'lon', 'ele', 'name')]

library(stringr)

Xpts$QuadName <- str_extract(Xpts$name, 'X[0-9]+')

plot_positions <- merge( Xplots, Xpts ) 

plot_positions <- subset( plot_positions , Treatment %in% c('Drought', 'Irrigation'))

plot_positions <- subset( plot_positions[ , c('QuadName', 'Treatment', 'lat', 'lon')], QuadName %in% c('X1', 'X7', 'X11', 'X15'))

plot_positions$plot <- paste( 'USSES', gsub( plot_positions$QuadName, pattern = 'X', replacement = ''), c(2, 12, 16, 16, 8), 'C' , sep = '_' )  

write.csv( file = 'data/control_plot_positions.csv', plot_positions[ , c('plot', 'lat', 'lon')], row.names = FALSE ) 
