library( ggplot2 )
library(dplyr)
library(tidyr)
library(stringr)
files <- dir('data/temp_data', 'recruitment.RDS', full.names = T)

dat <- lapply( files, readRDS )
rdat <- do.call(rbind, dat)

rdat <- subset(rdat, year > 2010 & Group == 'E1')


rratio <- 
  rdat %>% 
  gather( cover_spp, cover, starts_with('cov.')) %>% 
  mutate( cover_spp = str_extract(cover_spp , '[A-Z]{4}$')) %>% 
  filter( species == cover_spp ) 
  
ratplot <- 
  ggplot( rratio, aes( x = cover, y = Y)) + 
  geom_point() + 
  ylab( 'recruitment' ) + 
  xlab( 'parent species cover')  + 
  facet_grid(  Treatment ~ . ) + 
  theme_bw()

p_out <- rratio %>% 
  group_by( species ) %>% 
  do(p = ratplot %+% . + ggtitle(.$species))

pdf( 'figures/parent_species_cover_to_recruitment.pdf', height = 5, width = 5 ) 
print( p_out$p ) 
dev.off()

rrat_summary <- 
  rratio %>% 
  group_by( species, Treatment, Period ) %>% 
  filter( cover > 0 ) %>% 
  mutate( rratio = Y/cover ) %>% 
  summarise(rratio =  mean(rratio))

ggplot ( rrat_summary , aes ( x = Treatment, y = rratio ) ) + geom_point() + facet_grid ( species ~ . ) 


aggregate( data = rdat, Y ~ quad + Treatment + species , 'sum')

aggregate( data = rdat, Y ~ Treatment + species, 'mean')


base_plot <- 
  ggplot ( rdat , aes( x = Treatment, y = Y) ) + 
  geom_point(position = position_jitter(width = 0.2,height = 0), alpha = 0.5) + 
  ylab( 'Avg recruitment per plot in big exclosure 2011 to 2016') + 
  theme_bw()

pr <- 
  rdat %>% 
  group_by( species, Treatment, quad ) %>% 
  summarise( Y = mean(Y) ) %>% 
  group_by( species ) %>% 
  do(p =  base_plot %+% . + ggtitle( .$species))



rcov <- 
  rdat %>% 
  dplyr::select(-species, -Y )%>% 
  distinct() %>% 
  gather( species, cover, starts_with('cov')) %>% 
  group_by(species, Treatment, quad ) %>% 
  summarise( cover = mean(cover )) 

library(vegan)

mydat <- rcov %>% spread(species, cover ) %>% ungroup()
mds <- metaMDS(mydat[ , grep('cov', names(mydat))])
plot(mds, type = 't')
myfit <- envfit(mds, env = mydat[, 1] )
plot(myfit)

base_plot <- 
  ggplot ( rcov, aes( x = Treatment, y = cover) ) + 
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2,height = 0), alpha = 0.4) + 
  ylab( 'Average cover per quad 2011 to 2015') + 
  theme_bw()

pc <- rcov %>% group_by(species ) %>% do( p = base_plot %+% . + ggtitle(.$species))

pdf( 'figures/average_cover_in_experiment.pdf', width = 5, height = 5 )
print( 
  pc$p
)
dev.off()


library(gridExtra)

pdf( 'figures/total_experimental_recruitment.pdf', width = 6, height = 5 ) 

for( i in 1:length(pr$p)){ 
  print( 
    grid.arrange(pr$p[[i]], pc$p[[i]] , ncol = 2, nrow = 1) 
  )
}

dev.off()




rdat <- 
  rdat %>% 
  dplyr::select( year, quad, Treatment, species, Y) %>% 
  group_by( quad, species, Treatment ) %>% 
  summarise( Y = mean(Y) ) %>% 
  spread(species, Y)

pairs(rdat[, c('ARTR', 'HECO', 'POSE', 'PSSP')])

rdat <- rdat[ complete.cases(rdat) , ] 

cor(rdat[, c('ARTR', 'HECO', 'POSE', 'PSSP')])
mds <- metaMDS( comm = rdat[ ,c('ARTR', 'HECO', 'POSE', 'PSSP') ] )

plot(mds, type = 't')
myfit <- envfit(mds, env = mydat[, 1] )
plot(myfit)



