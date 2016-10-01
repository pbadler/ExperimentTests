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

summary(m)

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





