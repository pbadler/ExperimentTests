library(stringr)
library(gridExtra)
library(vegan)
library(dplyr ) 
library(tidyr)

quad_info <- read.csv('~/driversdata/data/idaho_modern/quad_info.csv')

get_cover_data <- 
  function(path = '~/driversdata/data/idaho/speciesData/') { 
  cover_files <- dir( path , 'quadratCover.csv', recursive = T, full.names = T)
  d <- lapply( cover_files, read.csv)
  species <- str_extract( cover_files, '[A-Z]{4}')
  d <-  mapply(x = d, y = species, function(x,y) {x$species <- y ; x }, SIMPLIFY = FALSE ) 
  d <- do.call(rbind, d)
  return( d )
}

historical_data <- get_cover_data()

modern_data <- get_cover_data('~/driversdata/data/idaho_modern/speciesData/')

historical_data$year <- historical_data$year + 1900

historical_data$Treatment <- 'Control'
historical_data$Period <- 'Historical'

modern_data <- merge(modern_data, quad_info[, c('quad', 'Treatment')], by = 'quad'  )
modern_data$Period <- 'Modern'

all_data <- rbind(historical_data, modern_data)

all_data <-
  all_data %>% 
  filter(! Treatment %in% c('No_grass', 'No_shrub')) %>% 
  group_by( Period) %>% 
  mutate( short_year = year - floor(median(year)))

cover_plot <- 
  ggplot (all_data, aes( x = year, y = log( totCover ), color = species, shape = Treatment, linetype = Treatment) ) + 
  geom_point(alpha = 0.1) + 
  stat_summary(fun.y = 'mean', geom = 'line', size = 1) +
  scale_y_continuous('log total cover', breaks = seq(-2, 8, by = 2), limits = c(-2, 9)) 

p1 <- cover_plot %+% subset(all_data, Period == 'Historical') + 
  scale_color_discrete(guide = FALSE) + 
  scale_linetype_discrete(guide = FALSE)+ 
  scale_shape_discrete(guide = FALSE) 

p2 <- cover_plot %+% subset(all_data, Period == 'Modern') + 
  xlim( 1998, 2025 ) + 
  theme(axis.title.y = element_blank())


print( grid.arrange(p1, p2 , ncol = 2, widths = c(0.44, 0.56) ) ) 

all_data_wide <- all_data %>% spread( species, totCover, fill = 0)
library(vegan)

cover_nmds <- metaMDS( all_data_wide[ ,c('ARTR', 'HECO', 'POSE', 'PSSP') ] ) 

fit <- envfit(cover_nmds, all_data_wide [ , c('Period', 'Treatment') ])
fit
fit
points( cover_nmds$points[all_data_wide$Period == 'Historical', ], col = 'red')
points( cover_nmds$points[all_data_wide$Period == 'Modern', ], col = 'black')
fit$factors
fit$vectors
