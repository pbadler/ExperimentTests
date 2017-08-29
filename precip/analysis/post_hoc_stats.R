rm(list = ls())
library(dplyr)
library(tidyr)

test <- read.csv('output/overall_lppd.csv')

files <- dir('data/temp_data/', pattern = 'survival.RDS', full.names = T)
all_dat <- lapply( files, readRDS)

survival <- do.call(rbind, all_dat)

nplants_survival <- survival %>% group_by(species, Period) %>% summarize( n = n())

nplants_survival
write.csv(nplants_survival, '~/Desktop/sample_size.csv')

library(ggplot2)

head(nplants_survival)
ggplot(nplants_survival, aes( x = species, y = n, fill = Period)) + geom_bar(stat = 'identity', position = 'dodge')

nplants_survival$vital_rate <- 'survival'

test <- test %>% spread(model, lppd) %>% mutate( diff = climate - year )

plot(test$diff)
nplants_survival
merge(test, nplants_survival %>% filter(Period == 'Historical'))

plot(data = merge(test, nplants_survival %>% filter(Period == 'Historical')), diff ~ n)

rm(list = ls())

read.cv('')

