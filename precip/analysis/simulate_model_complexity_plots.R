#
# simulate effect of model complexity on prediction quality 
#
#

# by model ---------------------------------------------------------------------------

x <- sort( rep(1:5, 10))
beta <- sort( rep( 10*c(1, 2, 3, 4, 5), 10 )) 
k <- sort( rep (  c(1, 10, 2, 12, 18), 10) ) # number of parameters 

y <- 100 + seq(-50, 50, length.out = 10) + x*beta

df <- data.frame(model = factor( x ) , lppd = y , k = k, model_label = factor(x, labels = c('1 - null', '2 - climate only', '3 - single species', '4 - single species + climate', '5 - multi-species + climate')))

model_types <- 
  ggplot( df, aes( x = model, y = lppd, fill = model_label ) ) + 
  geom_boxplot() + 
  xlab( 'Model' ) + 
  ylab( 'LPPD') + 
  scale_fill_discrete(name = '') + 
  theme( axis.text.y = element_blank(), 
         legend.text = element_text ( size = 16), 
         legend.title = element_text( size = 20), 
         axis.title = element_text(size = 20), 
         axis.text = element_text( size = 16), 
         legend.position = 'bottom', 
         legend.direction = 'vertical', 
         legend.title = element_blank()) 

png( 'figures/example_model_type_accuracy.png', 5, 5 , units= 'in', res = 300 )
print( model_types ) 
dev.off()

# by number of paramaters ------------------------------------------------------------

model_complexity <- 
  ggplot( df, aes( x = k , y = lppd) ) + 
  geom_point()  + 
  ylab( 'LPPD' ) + 
  xlab ( 'Number of parameters') + 
  theme( axis.text.y = element_blank(), 
         legend.text = element_text ( size = 16), 
         legend.title = element_text( size = 20), 
         axis.title = element_text(size = 20), 
         axis.text = element_text( size = 16)) 

png( 'figures/example_model_complexity_accuracy.png', 3.6, 4 , units= 'in', res = 300 )
print( model_complexity ) 
dev.off()
