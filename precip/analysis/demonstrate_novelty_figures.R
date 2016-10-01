#
# generate novel conditions figure for presentation 
#

library(ggplot2)


MM <- function( S, Vmax, K ) { Vmax*S / (K + S ) }
LL <- function( X, X0, L, k ){ L/(1 + exp(-k*(X-X0)))}


lm_right <- function(formula,data,...){
  # extend geom smooth past range of x values 
  mod <- lm(formula,data)
  class(mod) <- c('lm_right',class(mod))
  mod
}

predictdf.lm_left <- 
  # extend geom smooth past range of x values 
  function(model, xseq, se, level){
    init_range = range(model$model$x)
    ## here the main code: truncate to x values at the left
    xseq <- xseq[xseq <=init_range[2]]
    ggplot2:::predictdf.default(model, xseq[-length(xseq)], se, level)
  }

x <- seq(-20, 0, length.out = 100)
mu <- - LL(x, -15, 20, 1)
reps <- 2
y = rnorm( length(mu)*reps, rep(mu, 10), 3)

df <- data.frame( x = rep(x, reps),  y = y)

df$period <- NA
df$period[ df$x > -11] <- 'new' 
df$period[ df$x < -11] <- 'old'

base_plot <- 
  ggplot( df, aes( x = x, y = y ) ) + 
  geom_point( color = NA) + 
  xlab( expression( Seasonal ~ temperature ~ '('* degree * C * ')')) + 
  ylab( 'Ecological response') + 
  theme(axis.title.x =element_text( size = 20), 
        axis.title.y =element_text( size = 20), 
        axis.text.x = element_text( size = 15), 
        axis.text.y = element_blank()) + 
  xlim( c(-22, 0)) + 
  ylim( c(-40, 10))

historical_data <- 
  base_plot + 
  geom_point( data = subset(df, period == 'old') , aes( x = x, y = y )) 

historical_fit <- 
  base_plot + 
  geom_point( data = subset(df, period == 'old') , aes( x = x, y = y )) + 
  geom_smooth(  data = subset(df, period == 'old') , aes( x = x, y = y ), method = 'lm' ,se = FALSE )

predicted <- 
  base_plot + 
  geom_point( data = subset(df, period == 'old') , aes( x = x, y = y )) + 
  geom_smooth(  data = subset(df, period == 'old') , aes( x = x, y = y ), method = 'lm_right', fullrange = TRUE , se = FALSE) 

future <- 
  predicted + 
  geom_point(data = subset( df, period == 'new'), aes( x = x, y = y ) , color = 'red') 

png( 'figures/novelty_base_plot.png', width = 6, height = 4.8, units = 'in', res = 300)
print( base_plot )  
dev.off()

png( 'figures/novelty_historical_data.png', 6, 4.8, 'in', res = 300)
print ( historical_data ) 
dev.off()

png( 'figures/novelty_historical_fit.png', 6, 4.8, 'in', res = 300)
print ( historical_fit ) 
dev.off()

png( 'figures/novelty_predicted.png', 6, 4.8, 'in', res = 300)
print ( predicted ) 
dev.off()

png( 'figures/novelty_future.png', 6, 4.8, 'in', res = 300)
print ( future ) 
dev.off()

