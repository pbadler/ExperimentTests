rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)

climate <- readRDS('data/temp_data/all_clim_covs.RDS')

#lppd_scores <- read.csv('output/lppd_scores.csv')
error <- read.csv('output/climate_effect_prediction_error.csv')
datlist <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')

# base plot ------------------------------------------------------------- 

cover <- unique( data.frame( year = datlist$PSSP$year2, Treatment = datlist$PSSP$treat2, Group = datlist$PSSP$gid2, quad = datlist$PSSP$quad2, datlist$PSSP$parents12))

cover$Period <- NA
cover$Period[ cover$year < 2006] <- 'Historical'
cover$Period[ cover$year > 2006] <- 'Modern'

cover$Treatment <- factor(cover$Treatment, labels = c('Control', 'Drought', 'Irrigation'))

ggplot(cover, aes( x = year, y = cov.ARTR, group = quad)) + geom_point() + geom_line() + geom_smooth( aes( group = Period)) + facet_grid( Group ~ . )
ggplot(cover, aes( x = year, y = cov.HECO, group = quad)) + geom_point() + geom_line() + geom_smooth( aes( group = Period)) 
ggplot(cover, aes( x = year, y = cov.POSE, group = quad)) + geom_point() + geom_line() + geom_smooth( aes( group = Period))
ggplot(cover, aes( x = year, y = cov.PSSP, group = quad)) + geom_point() + geom_line() + geom_smooth( aes( group = Period))

cover <- 
  cover %>% 
  gather( species, cover.t0, starts_with('cov')) 

cover.t1 <- cover %>%
  mutate( year = year - 1 ) %>% 
  rename( cover.t1 = cover.t0 )

cover <- merge( cover, cover.t1) %>% arrange( species, quad, year ) 

climate  <- 
  climate %>% 
  mutate( P.a.t3 = P.a.1 + P.a.0 + P.a.l , T.g.t3 = (T.sp.1 + T.sp.0 + T.sp.l + T.su.1 + T.su.0 + T.su.l + T.f.1 + T.f.0 + T.f.l)/9 )

climate[ , grep( '(^P\\.)|(^T\\.)|(^VWC\\.)', names( climate ) ) ] <- scale( climate[ , grep( '(^P\\.)|(^T\\.)|(^VWC\\.)', names( climate ) ) ] )

cover <- merge(cover, climate)

library(lme4)

cover$logcover.t1 <- log(cover$cover.t1)
cover$logcover.t0 <- log(cover$cover.t0)
cover$lc <- cover$logcover.t1 - cover$logcover.t0


ggplot(subset( cover, species == 'cov.ARTR'), aes( x = year, y = lc)) + geom_point() + geom_smooth( aes( group = Period)) + facet_grid( Group ~ . )
ggplot(subset( cover, species == 'cov.HECO'), aes( x = year, y = lc)) + geom_point() + geom_smooth( aes( group = Period)) + facet_grid( Group ~ . )
ggplot(subset( cover, species == 'cov.POSE'), aes( x = year, y = lc)) + geom_point() + geom_smooth( aes( group = Period)) + facet_grid( Group ~ . )
ggplot(subset( cover, species == 'cov.PSSP'), aes( x = year, y = lc)) + geom_point() + geom_smooth( aes( group = Period)) + facet_grid( Group ~ . )



m1 <- lmer(data = subset( cover, species == 'cov.PSSP' & is.finite( logcover.t0) & is.finite(logcover.t1)), logcover.t1 ~ logcover.t0 + factor(Group) + Period + factor(Treatment) + (1|year))
summary(m1)

m2 <- update(m1, . ~ . + P.a.t3*T.g.t3)
summary( m2 ) 




plot(cover$year, cover$cov.ARTR)



library(vegan)
temp_pca <- metaMDS(cvs_in)
dists <- apply( rbind( cvs_in, cvs_out), 1, function( x ) dist(rbind( rep(0, length(x)), x) ))

temp_pca_both <- metaMDS(rbind(cvs_in, cvs_out), trymax = 100)

plot_df <- 
  data.frame( temp_pca_both$scores,
                       Period = cover$Period,                        
                       Treatment = cover$Treatment, 
                       year = cover$year, 
                       dists = dists)

plot_df$Treatment <- factor( plot_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))

biplot(temp_pca_both)

plot_df <- 
  plot_df %>% 
  mutate(center1 = mean(Comp.1[Period == 'Historical']), center2 = mean(Comp.2[Period == 'Historical']))

df_arrows  <- data.frame( temp_pca_both$loadings[ , c(1:2)], x = temp_pca_both$center[1], y = temp_pca_both$center[2])

names( df_arrows ) [ 1:2 ] <- c('xend', 'yend')

df_arrows$xend <- df_arrows$xend*max(plot_df$Comp.1)
df_arrows$yend <- df_arrows$yend*max(plot_df$Comp.2)

gg_PCA <- 
  ggplot( plot_df, aes( x = Comp.1, y = Comp.2, color = Period )) + 
  geom_point(size = 2) + 
  geom_segment(data = df_arrows, aes (xend = xend, yend = yend, x = x, y = y), color = 1, arrow = arrow()) + 
  geom_text( data = df_arrows , aes( x = xend*1.5, y = yend*1.5, label = row.names(df_arrows)), color = 1) 


gg_PCA <- make_gg_PCA('ARTR', datlist)  
png(filename = 'figures/climate_pca_comparison.png', height = 6, width = 6, units = 'in', res = 300)
print( gg_PCA ) 
dev.off()

# plot novelty of years -------------------------------------------------

pca_dat <- gg_PCA$data

clim <- unique( data.frame( year = datlist$ARTR$year2 , Treatment = datlist$ARTR$treat2, datlist$ARTR$C2) ) 
clim$Treatment <- factor(clim$Treatment, labels = c('Control', 'Drought', 'Irrigation'))
pca_dat <- merge( pca_dat, clim, by = c('year', 'Treatment') )

my_colors <- c('#1b9e77', '#d95f02', '#7570b3')

ggplot(subset( pca_dat, year > 2005), aes(x = year, y = dists , color = Treatment )) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = my_colors)

ggplot(subset( pca_dat, year > 2005), aes(x = year, y = T.sp.1)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = my_colors)

ggplot(subset( pca_dat, year > 2005), aes(x = year, y = VWC.sp.1, color = Treatment)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = my_colors)


# plot lppd by distance -------------------------------------------------
gg_PCA <- make_gg_PCA('ARTR', datlist)

plot_df <- gg_PCA$data

MEA <- 
  error %>% 
  filter( type == 'error') %>%
  group_by(vital_rate, species, Treatment, year, type  ) %>% 
  summarise( MAE = mean(abs(val)) , ME = mean(val), sd_e = sd(val))

error_dist <- merge( plot_df, MEA , by = c('Treatment', 'year')) 

ggplot( error_dist, aes( x = dists, y = MAE)) + 
  geom_point() + 
  geom_smooth( method = 'lm', group = 1, se = F)

ggplot( error_dist, aes( x = dists, y = MAE)) + geom_point() + 
  facet_grid( ~ species ) + 
  geom_smooth( method = 'lm', group = 1, se = F)

ggplot( error_dist, aes( x = dists, y = abs( ME ) )) + geom_point() + 
  facet_grid( ~ species ) + 
  geom_smooth( method = 'lm', group = 1, se = F)

ggplot( error_dist, aes( x = dists, y = sd_e )) + geom_point() + 
  facet_grid( ~ species ) + 
  geom_smooth( method = 'lm', group = 1, se = F)



lppd_scaled_scores <- 
  lppd_scores %>% 
  filter( Period == 'Modern') %>% 
  group_by(vital_rate, species ) %>% 
  mutate( lppd_scaled = scale(lppd))

lppd_scaled_scores$year_id <- as.numeric ( factor( lppd_scaled_scores$year ) )

plot_df$treatment <- factor( plot_df$treatment, labels = c('Control', 'Drought', 'Irrigation'))

test <- subset( plot_df, Period =='Modern')


lppd_scaled_scores$lppd_scaled <- as.numeric(lppd_scaled_scores$lppd_scaled)

lppd_scaled_scores <- merge( lppd_scaled_scores, test, by = c('treatment', 'year_id'))

plot_dists <- 
  lppd_scaled_scores %>% 
  group_by(dists) %>% 
  summarise(mlppd = mean(lppd_scaled), ucl = quantile(lppd_scaled, 0.9), lcl = quantile( lppd_scaled, 0.1))

gg_dists <- 
  ggplot(lppd_scaled_scores, aes( x = dists, y = lppd_scaled  )) + 
  geom_errorbar(data = plot_dists, aes(x = dists, y = mlppd, ymin = lcl, ymax = ucl), alpha = 0.3)  + 
  geom_point( data = plot_dists, aes( x = dists, y = mlppd) )  


png(filename = 'figures/climate_dists.png', heigh = 4, width = 4, units = 'in', res = 300)
print( gg_dists + ylab( 'relative predictive score') + xlab( '"novelty"') ) 
dev.off()

