#
# Simulate predictive error
# generate plots showing possible responses 
#
#
rm(list = ls())
library(ggplot2)

x <- rnorm( 10000)
x <- c(x, rnorm(10000, -2 ))

fit <- data.frame(x = x, type = sort( rep(c('Historical', 'Modern'), 0.5*length(x) )))

theme_set ( theme_get( ) + theme( strip.text = element_text(size = 12), plot.title = element_text(size = 14 ))  ) 

ex1 <- 
  ggplot ( data = fit, aes(x = x , fill = type) ) + 
  geom_density(color = NA, alpha = 0.5) + 
  geom_vline(aes(xintercept = 0)) + 
  facet_wrap(~type, ncol = 1) + 
  ylab( "Density") + 
  xlab( "Size Predicted - Observed") + 
  scale_fill_discrete(guide = FALSE) + 
  xlim ( -10, 10)

x <- x[1:10000] 
x <- c(x, rnorm(10000, 0, 3))

fit2 <- data.frame(x = x, type = sort( rep(c('Historical', 'Modern'), 0.5*length(x) )))

ex2 <- ex1 %+% fit2

png( 'figures/example_growth_error_1.png', 3.5, 3.5, units = 'in', res = 300)
print( ex1  + ggtitle( "Growth") ) 
dev.off()

png( 'figures/example_growth_error_2.png', 3.5, 3.5, units = 'in', res = 300)
print( ex2  + ggtitle( "Growth") ) 
dev.off()

# -- for survival predictions -------------------------------------------------------------------------- # 
rm(list = ls())

library(ggplot2)
library(ROCR)
library(caret)
library(pROC)

# functions ---------------------------------------- # 

inv_logit <- function( x ) exp(x)/(1 + exp(x))

# parameters --------------------------------------- # 

x <- runif( 1000, -10, 10)
beta1 <- 1
beta2 <- 0.5
sig_b1 <- c(0.5, 0.1)
sig_b2 <- c(0.5, 0.1)

# model -------------------------------------------- # 
df <- data.frame(training = 0:1, x = x, beta1 = beta1, beta2 = beta2, sig_b1 = sig_b1, sig_b2 = sig_b2 )
df$mu <- rnorm(length(x), df$beta1, df$sig_b1) + rnorm(length(x), df$beta2, df$sig_b2)*df$x  

# random ------------------------------------------- # 

df$p <- inv_logit ( df$mu )
df$y <- rbinom(length(x), 1, df$p)

# organize ----------------------------------------- #
training <- subset(df, training == 1)
hold_out <- subset(df, training == 0)

# plot --------------------------------------------- # 
print( 
  ggplot( df, aes ( x = x, y = y , color = factor( training) )) + 
    geom_count(position = position_jitter(height = 0.1)) + 
    geom_smooth(method = 'glm', method.args = list( family = 'binomial')) )

# fit new model ------------------------------------ # 

m <- glm(data = training, y ~ x, family = 'binomial' ) 

# predict -------------------------------------------# 
prob_out <- predict( m , newdata = hold_out, type = 'response')
prob_train <- predict( m , newdata = training, type = 'response')

# evaluate with ROC and AUROC ---------------------- # 
pred_train <- prediction(prob_train, training$y )
pred_out <- prediction( prob_out, hold_out$y )

perf_train <- performance(pred_train, measure = "tpr", x.measure = "fpr")
plot(perf_train)

perf_out <- performance( pred_out, measure = 'tpr', x.measure = 'fpr')
plot(perf_out)

auc_train <- performance(pred_train, measure = "auc")
auc_train <- auc_train@y.values[[1]]
auc_train

auc_out <- performance(pred_out, measure = "auc")
auc_out <- auc_out@y.values[[1]]
auc_out

length(perf_out@y.values[[1]])

plot_roc <-data.frame(x <- c(perf_train@x.values[[1]], perf_out@x.values[[1]]),  y <- c(perf_train@y.values[[1]], perf_out@y.values[[1]]) , type <- sort( rep( c('Historical', 'Modern'), length(perf_out@y.values[[1]]) )))

png ( 'figures/example_survival_error.png', 3.5, 3.5, units = 'in', res = 300 )
print ( 
  ggplot( data = plot_roc, aes( x = x, y = y, color = type )) + 
    geom_step(size = 1.5) + 
    xlab( 'False Positive Rate') + 
    ylab( 'False Negative Rate') +
    theme(legend.title = element_blank(), legend.text = element_text(size  = 12) , legend.position = c(0.75, 0.25) ) +
  ggtitle( "Survival")  
)
dev.off()

# -- for recruitment predictions ----------------------------------------------------------------------------------- # 

x <- rpois( 10000, lambda = 20) - 20 

x <- c(x, rpois(10000, lambda = 10 ) - 20 ) 

fit <- data.frame(x = x, type = sort( rep(c('Historical', 'Modern'), 0.5*length(x) )))

ex1 <- 
  ggplot ( data = fit, aes(x = x , fill = type) ) + 
  geom_histogram(alpha = 0.5) + 
  geom_vline(aes(xintercept = 0)) + 
  facet_wrap(~type, ncol = 1) + 
  xlab( "No. Recruits, Predicted - Observed") + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) + 
  scale_fill_discrete(guide = FALSE) 

png( 'figures/example_recruitment_error.png', 3.11, 3.5, units = 'in', res = 300)
print( ex1 + ggtitle( "Recruitment") ) 
dev.off()

