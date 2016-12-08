climate_covariates <- readRDS('data/temp_data/all_clim_covs.RDS')

climate_covariates <- climate_covariates[ complete.cases(climate_covariates), ]

climate_covariates <- subset( climate_covariates , !(Treatment != 'Control' & year < 2011 ))

h <- subset( climate_covariates, Period == 'Historical' & Treatment == 'Control')
m <- subset( climate_covariates, Period == 'Modern')

cn <- grep( '^[VTP].*\\..*1$', names( climate_covariates))
cn <- cn[-2]

hscaled <- scale(h[,cn])
h[,cn] <- hscaled
m[, cn] <- scale(m[, cn], attr(hscaled, 'scaled:center'), attr(hscaled, 'scaled:scale'))

hcor <- cor(h[, cn])
mcor <- cor(h[, cn])

hcov <- cov(h[, cn])
mcov <- cov(h[, cn])

colMeans(h[,cn])
apply( h[,cn], 2, sd)

library(mvtnorm)

x <- seq(-2, 2, length.out = 1000)
plot(x,  dnorm(x, mean = 0, sd = 1))
h[1,cn]

m$ll <- dmvnorm(m[,cn], sigma = hcov, log = T)
h$ll <- dmvnorm(h[, cn], sigma = hcov, log = T)

plot(m$year, m$ll, col = 1)
points(m$year[m$Treatment == 'Drought' ], m$ll[m$Treatment == 'Drought'], col = 'red')
points(m$year[m$Treatment == 'Irrigation' ], m$ll[m$Treatment == 'Irrigation'], col = 'blue')

plot(h$year, h$ll, col = 1)
points(h$year[h$Treatment == 'Drought' ], h$ll[h$Treatment == 'Drought'], col = 'red')
points(h$year[h$Treatment == 'Irrigation' ], h$ll[h$Treatment == 'Irrigation'], col = 'blue')
