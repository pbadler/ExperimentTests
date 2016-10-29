rm(list = ls())

library(lme4)

true_cover <- read.csv('~/driversdata/data/idaho_modern/speciesData/ARTR/quadratCover.csv')
quad_info  <- read.csv('~/driversdata/data/idaho_modern/quad_info.csv')
true_cover <- merge(true_cover, quad_info)
true_cover$cover <- 100*(true_cover$totCover/10000)
true_cover$year <- true_cover$year 

true_cover_per_treatment <- 
  true_cover %>% 
  filter( !Treatment %in% c('No_shrub','No_grass')) %>% 
  group_by(Treatment, year ) %>% 
  mutate(cover = mean(cover))

# model growth 

gdf <- readRDS('data/temp_data/ARTR_scaled_growth_dataframe.RDS')

gdf <- subset(gdf, Period == 'Modern')

C <- names(gdf)[grep('^[PT]\\.', names(gdf))]

base_f <- as.formula( Y ~ (X|yid) + Group + W.ARTR )

f <- as.formula( paste( c( expression(Y ~ X + (X|yid) + Group + W.ARTR), paste(C,collapse = ' + ')), collapse = ' + '  )  ) 

gm1 <- lmer(data = gdf, base_f )
summary(gm1)

gdf$yhat <- predict(gm1)
plot(gdf$yhat, gdf$Y)
abline(0, 1)

# model survival 

sdf <- readRDS('data/temp_data/ARTR_scaled_survival_dataframe.RDS')

sdf <- sdf[ sdf$yid %in% gdf$yid , ]

sdf <- subset(sdf, Period = 'Modern')

C <- names(gdf)[grep('^[PT]\\.', names(gdf))]

f <- as.formula( paste( c( expression(Y ~ X + (X|yid) + Group + W.ARTR), paste(C,collapse = ' + ')), collapse = ' + '  )  ) 

sm1 <- glmer(data = sdf, base_f , family = 'binomial')
summary(sm1)

sdf$yhat <- predict(sm1, type = 'response' )
plot(sdf$yhat, sdf$Y)
abline(0, 1)

# 
cover_df <- sdf 
cover_df$size_predicted <- predict( gm1, sdf )
cover_df$survival_predicted <- predict(sm1, sdf, type = 'response')

total_cover_pred <- 
  cover_df %>% 
  mutate( year = year + 1 ) %>%                                         # predictions are for the next year 
  mutate( size_adjusted = exp(size_predicted)*survival_predicted) %>% 
  group_by(year, quad) %>% 
  summarise( cover = 100*(sum(size_adjusted)/10000))

total_cover_pred <- merge( total_cover_pred, quad_info)

cover_per_treatment_pred <-  
  total_cover_pred %>% 
  group_by(year, Treatment) %>% 
  summarise( cover = mean(cover ) ) 

ggplot( cover_per_treatment_pred, aes(x = year, y = cover, color = Treatment )) + 
  geom_line(linetype = 2 ) + 
  geom_line(data = true_cover_per_treatment, aes( x = year, y = cover, color = Treatment )) + 
  facet_wrap(~Treatment)

