rm(list = ls())
library(lme4)

mydf <- readRDS('data/temp_data/ARTR_scaled_growth_dataframe.RDS')

# these are outputs from Chengjin's data processing script 
# ARTR_survival_gaussW.r in the Nature Comm supplemental information 
ccdf <- read.csv('data/temp_data/cc_ARTR_surv.csv')
ccW <- read.csv('data/temp_data/cc_Ws.csv')
mydf <- subset(mydf, Period == 'Historical')

table(ccdf$survives, ccdf$year)
table(mydf$survives, mydf$year)

mean(ccdf$logarea)
mean(mydf$logarea.t0)

ccdf$TmeanSpr1_scaled <- scale(ccdf$TmeanSpr1)
ccdf$pptLag_scaled <- scale(ccdf$pptLag)
ccdf$year <- ccdf$year + 1900 

ccdf_temp <- unique(ccdf[ , c('year', 'TmeanSpr1_scaled')])
my_temp <- unique(mydf[, c('year', 'T.sp.0')])

ccdf_ppt <- unique(ccdf[, c('year', 'pptLag_scaled')])
my_ppt <- unique(mydf[ , c('year', 'P.f.w.sp.l')])

plot( my_temp, type = 'l')
points( ccdf_temp, type = 'l', col = 'red')

plot(my_ppt, type = 'l')
points(ccdf_ppt, type = 'l', col = 'red')

mymod <- glm(data = mydf, survives ~ Group+ logarea.t0*P.f.w.sp.l + W.ARTR, family = 'binomial')
ccmod <- glm(data = ccdf, survives ~ Group + logarea*pptLag_scaled + ccW[,1 ], family = 'binomial')

summary(mymod)
summary(ccmod)

head(ccW)
head(mydf$W.ARTR)

nrow(ccW)
nrow(mydf)

summary( mydf$P.f.w.sp.l)
summary( mydf$`P.f.w.sp.l:logarea.t0`)

mymod <- glm(data = mydf, survives ~ Group+ logarea.t0 + P.f.w.sp.l + `P.f.w.sp.l:logarea.t0` + W.ARTR, family = 'binomial')
summary(mymod)

mymodmer <- glmer( data = mydf, survives ~ Group + logarea.t0 + P.f.w.sp.l + `P.f.w.sp.l:logarea.t0` + W.ARTR + (logarea.t0|year), family = 'binomial')
summary(mymodmer)

ccdf$`pptLag:logarea_scaled` <- scale( ccdf$pptLag*ccdf$logarea ) 

ccmodmer <- glmer( data = ccdf, survives ~ Group + logarea + pptLag_scaled + `pptLag:logarea_scaled` + ccW[, 1] + (logarea|year), family = 'binomial')
summary(ccmodmer)

# chengjin and my growth data 

ccdf <- read.csv('data/temp_data/cc_ARTR_growth.csv')
ccdfW <- read.csv('data/temp_data/cc_ARTR_Ws_growth.csv')

mydf <- readRDS('data/temp_data/ARTR_scaled_growth_dataframe.RDS')

mydf <- subset(mydf, Period == 'Historical')
nrow(ccdf)
nrow(mydf)

ccdf$year <- ccdf$year + 1900 

all_df <- merge( ccdf, mydf, by = c('trackID', 'year', 'quad', 'Group'), all.x = TRUE)

all_df[which(is.na(all_df$logarea.t0.y)), ]

plot(all_df$logarea.t1.x, all_df$Y)
plot(all_df$TmeanSpr2, all_df$T.sp.1)
plot(all_df$pptLag, all_df$P.a.l)

ccmod <- lmer(data = all_df, logarea.t1.x ~ logarea.t0.x + (1|Group) + (logarea.t0.x|year) + W.ARTR + W.POSE + W.HECO + W.PSSP + TmeanSpr2 + TmeanSpr1 + ppt1 + ppt2 + pptLag)
mymod <- lmer(data = all_df, logarea.t1.y ~ logarea.t0.y + (1|gid) + (logarea.t0.y|year) + W.ARTR + W.POSE + W.HECO + W.PSSP + T.sp.0 + T.sp.1 + P.f.w.sp.0 + P.f.w.sp.1 + P.a.l)

summary(ccmod) 
summary(mymod)


