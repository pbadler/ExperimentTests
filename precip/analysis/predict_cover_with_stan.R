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

base_f <- as.formula( Y ~ (X|yid) + Group + W.ARTR + T.sp.0 + T.sp.1)

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

# model survival stan 

N <- nrow(sdf)
Y <- sdf$Y
X <- sdf$X
yid <- sdf$yid - 22
nyrs <- nlevels(factor(yid))
G  <- nlevels(sdf$gid)
gm <- model.matrix.lm(~sdf$gid)
gid <- as.numeric(sdf$gid)
W <- sdf$W.ARTR
Wcovs <- 1 
sdatalist <- list( N = N, Y = Y, X = X, yid = yid , nyrs =nyrs, W = W, Wcovs = Wcovs, gid = gid, gm = gm)

sfit <- stan('analysis/survival/model_survival_1.stan', data = sdatalist, pars = c('mu', 'sig_a', 'sig_b1', 'w', 'bg'), chains = 1, iter = 1000)

sdf$s_mu <- summary(sfit, 'mu')$summary[, 1] 

# model growth stan 

N <- nrow(gdf)
Y <- gdf$Y
X <- gdf$X
yid <- gdf$yid - 22
nyrs <- nlevels(factor( yid ))
G  <- nlevels(gdf$gid)
gid <- as.numeric(gdf$gid)
gm <- model.matrix.lm(~gdf$gid)
W <- gdf$W.ARTR
Wcovs <- 1 

gdatalist <- list( N = N, Y = Y, X = X, yid = yid , nyrs =nyrs, W = W, Wcovs = Wcovs, gid = gid, gm = gm )
datalist3 <- sdatalist
names(datalist3) <- paste0( names( datalist3) , 3 )
gdatalist <- c(gdatalist, datalist3 )

gfit <- stan('analysis/growth/model_growth_1_cover_predict.stan', data = gdatalist, pars = c('mu_out', 'sigma', 'sig_a', 'sig_b1', 'w'), chains = 1, iter = 1000)

sdf$g_mu <- summary(gfit, 'mu_out')$summary[, 1]

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

# plot stan predictions 

total_cover_stan <- 
  cover_df %>% 
  mutate( year = year + 1 ) %>%
  mutate( size_adjusted_stan = exp(g_mu)*s_mu ) %>% 
  group_by( year, quad) %>% 
  summarise( cover = 100*(sum(size_adjusted_stan)/10000))

total_cover_pred_stan <- merge(total_cover_stan, quad_info)

cover_per_treatment_stan <- 
  total_cover_pred_stan %>% 
  group_by(year, Treatment) %>% 
  summarise( cover = mean(cover ))

ggplot( cover_per_treatment_stan, aes( x = year, y = cover, color = Treatment)) + 
  geom_line(linetype = 2) + 
  geom_line(data = true_cover_per_treatment, aes( x = year, y= cover, color = Treatment ) ) + 
  facet_wrap(~Treatment)

#
gdf %>% 
  group_by( Group , Treatment  ) %>% distinct()

levels(gdf$Group)

summary(gfit, 'bg')$summary
summary(gm1)

df2 <- readRDS('data/temp_data/ARTR_growth.RDS')

clim <- readRDS('data/temp_data/all_clim_covs.RDS')
clim$year <- clim$year - 1 

df2 <- merge(df2, clim , by = c('Treatment', 'year', 'Period' ))

df_old <- subset(df2, Period == 'Historical')
m2_old <- lmer(data = df_old, logarea.t1 ~ logarea.t0  + (logarea.t0|year) + T.sp.0 + T.sp.1 + T.sp.1:logarea.t0 + T.sp.0:logarea.t0)
summary(m2_old)

df_new <- subset(df2, Period == 'Modern')
m2_new <- lmer(data = df_new, logarea.t1 ~ logarea.t0 + (logarea.t0|year))
summary(m2_new)
fixef(m2_new)
fixef(m2_old)

df <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')[['ARTR']]
df$W <- df$W[,1]
gfit_hist <- stan('analysis/growth/model_growth_1.stan', data = df, pars = c('mu', 'bg') , chains = 1, iter = 1000 )

summary(gfit_hist, 'bg')$summary
summary(gm1)

df1 <- readRDS('data/temp_data/ARTR_scaled_growth_dataframe.RDS')
df2 <- readRDS('data/temp_data/ARTR_growth.RDS')
clim <- readRDS('data/temp_data/all_clim_covs.RDS')
clim$year <- clim$year 
df2 <- merge(df2, clim, by = c('Treatment', 'year'))

df_both <- left_join( df1, df2 , by = c('year', 'quad', 'Treatment', 'trackID', 'species', 'Group')) 
df_both <- subset(df_both, Period == 'Historical')

plot( df_both$T.sp.0.x, df_both$T.sp.0.y)
plot( df_both$`T.sp.0:logarea.t0`, df_both$T.sp.0.y*df_both$logarea.t0.y)

m1 <- lm(data = df_both, logarea.t1.x ~ logarea.t0.x + T.sp.0.x + `T.sp.1:logarea.t0` )
m2 <- lm(data = df_both, logarea.t1.x ~ logarea.t0.x + T.sp.0.y + T.sp.1.y:logarea.t0.y )

summary(m1)
summary(m2)

plot(df_both$year, df_both$T.sp.0.x, type = 'l', ylim = c(-3, 24))
points(df_both$year, df_both$T.sp.0.y, type ='l', col = 'red')
clim2 <- subset(clim, Period == 'Historical')
points(clim2$year, clim2$T.sp.0, type = 'l', col = 'blue')

