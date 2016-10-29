growth_df <- readRDS('data/temp_data/ARTR_growth.RDS')
clim <- readRDS('data/temp_data/all_clim_covs.RDS')

clim_vars <- c( 'P.a.l', 'P.f.w.sp.0', 'P.f.w.sp.1', 'P.su.0', 'P.su.1', 'T.sp.0', 'T.sp.1')
clim <- clim[ , c('year', 'Treatment', 'Period', clim_vars)]

df <- merge(growth_df, clim, by = c('year', 'Treatment', 'Period'))

ifx <- df[, clim_vars]*df[, 'logarea.t0']   ##### climate by size interactions
names(ifx ) <- paste0(clim_vars , 'x', 'logarea.t0')
df <- cbind(df, ifx)

covars <- grep( '^[PTW]\\.', names (df ) ) 
scaled <- scale(df[ , covars])
df[, covars] <- scaled 

#if ( vr == 'growth' ){ 
  df$X <- scale(df$logarea.t0)
  df$Y <- scale( df$logarea.t1 ) 
#  } else if ( vr == 'survival'){
#  df$X <- scale(df$logarea.t0)
#  df$Y <- df$survives
#} else if ( vr == 'recruitment'){ 
#  df$Y <- df$Y
#  }

m1 <- lmer( data = df, Y ~ X + (X|year) + T.sp.0 + T.sp.0xlogarea.t0 )

# unscaled 
growth_df <- readRDS('data/temp_data/ARTR_growth.RDS')
clim <- readRDS('data/temp_data/all_clim_covs.RDS')

clim_vars <- c( 'P.a.l', 'P.f.w.sp.0', 'P.f.w.sp.1', 'P.su.0', 'P.su.1', 'T.sp.0', 'T.sp.1')
clim <- clim[ , c('year', 'Treatment', 'Period', clim_vars)]

df <- merge(growth_df, clim, by = c('year', 'Treatment', 'Period'))

ifx <- df[, clim_vars]*df[, 'logarea.t0']   ##### climate by size interactions
names(ifx ) <- paste0(clim_vars , 'x', 'logarea.t0')
df <- cbind(df, ifx)
df$Y <- df$logarea.t1
df$X <- df$logarea.t0

m2 <- lmer( data = df, Y ~ X*Grazing + (X|year) + T.sp.0 + T.sp.0xlogarea.t0 )

summary(m1)
summary(m2)

plot( predict(m2), predict(m1) ) 

sum( df$W.ARTR ) 
hist(df$W.ARTR)
hist(scale(df$W.ARTR))
hist(scale(scale(df$W.ARTR)))

