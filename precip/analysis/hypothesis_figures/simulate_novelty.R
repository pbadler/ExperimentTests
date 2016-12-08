b0 = 0 
b1 = 1
X <- c(1:100)
Y <- X

for( i in X ) { 
  Y[i] <- b0 + b1*X[i] + rnorm( 1, 1, 10 )
}

m <- lm( Y ~ X )

mu <- predict(m)

plot( Y)
points( mu , type = 'l', col = 'red')


