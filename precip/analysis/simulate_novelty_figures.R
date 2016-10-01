#
# Simulate novelty figure 
#
#

rm(list = ls()) 
library(ggplot2)

theme_set( theme_get() + theme( strip.text = element_text(size = 20), 
       axis.title = element_text(size = 20), 
       axis.text = element_text( size = 16), 
       legend.text = element_text(size = 16), 
       legend.title = element_text(size = 20)))

my_colors <- c( 'gray', '#fc8d62', '#8da0cb' )

x <- sort( rep( c(1,2,3), 1000))

X <-model.matrix.default( ~ factor(x))
beta <- c(100, -10, -10 )

y <- rnorm(length(x),  X%*%beta)

df <- data.frame( x, y , treatment = factor( x, labels = c('control', 'drought', 'irrigation')), novelty = factor( x, labels = c('historical', 'novel', 'novel')))

treatment_novelty_plot <- 
  ggplot( df, aes( x = y, fill = treatment)) + 
  geom_density(color = NA, alpha = 0.7) + 
  scale_fill_manual(values = my_colors, guide = FALSE ) + 
  ylab ( '') + 
  xlab ( 'log pointwise predictive density (LPPD)') + 
  facet_wrap( ~ treatment, nrow = 1 ) + 
  coord_flip()


# ----plot by multivariate novelty -----------------------------------------------------------------------------# 

x <- rnorm( 400 ) 
y <- rnorm( 400 ) 

df <- data.frame(x, y, period = 'Historical')

x <- rnorm ( 100, 2) 
y <- rnorm ( 100, 2, 0.5)

dfmod <- data.frame(x, y, period = 'Modern')

df <- rbind( df, dfmod ) 
df_centroid <- data.frame( t(colMeans(df[, c('x', 'y')])) ) 

names( df_centroid ) <- c('center.x', 'center.y')

df <- cbind( df, as.data.frame( df_centroid ) )

arrow_df <- head( subset(df, period == 'Modern'), 10 ) 

pca_plot1 <- 
  ggplot( df, aes(x , y , color = period ) ) + 
  geom_point() + 
  xlab ('PC1') + 
  ylab ('PC2') + 
  scale_color_manual(values = c('darkgray', my_colors[2]) ) +
  theme( legend.title = element_blank(), legend.position = 'bottom', legend.direction = 'vertical')

pca_plot2 <- pca_plot1 + 
  geom_segment(data = arrow_df, aes (x = center.x, y = center.y, yend = y, xend = x ), color = 'black') 


# now plot simple predictive score vs. novelty 

set.seed(2)
d <- exp(seq(0.1, 4, length.out = 100))

mu <- 1/d
y <- rnorm( length(d), mu, 0.4/d )

df <- data.frame( d, y ) 

lppd_by_novelty_plot <- 
  ggplot( data = df , aes( x = d, y = y ) ) + 
  geom_point() + 
  ylab('LPPD') + 
  xlab('Novelty\n(distance from historical centroid)') + 
  theme(axis.text = element_blank() )

lppd_by_novelty_plot

# save plots ------------------------------------------------------------------------------------ # 

png( 'figures/treatment_novelty_plot.png', 6, 6, units = 'in', res = 300)
print(treatment_novelty_plot)
dev.off()

png( 'figures/novelty_pca_plot1.png', 4.5, 5, units = 'in', res = 300)
print( pca_plot1 ) 
dev.off()

png( 'figures/novelty_pca_plot2.png', 4.5, 5, units = 'in', res = 300)
print( pca_plot2 ) 
dev.off()

png('figures/lppd_by_novelty_plot.png', 4.5, 5, units = 'in', res = 300)
print( lppd_by_novelty_plot ) 
dev.off()
