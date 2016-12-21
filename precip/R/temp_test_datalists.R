rm(list = ls())
gtest <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[[1]]
stest <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[[1]]

gdf <- readRDS('data/temp_data/ARTR_growth_cleaned_dataframe.RDS')
sdf <- readRDS('data/temp_data/ARTR_survival_cleaned_dataframe.RDS')

testdf <- merge( sdf, gdf , by = c('trackID', 'year', 'quad', 'Group'))
plot( testdf[, c('X.x', 'X.y')] ) # matches 
dim(sdf)
dim(gdf)
dim(testdf)
head(testdf)

test1 <- data.frame( track = gtest$trackidhold, year = gtest$yearhold, quad = gtest$quadhold, Xhold  = gtest$Xhold, group = gtest$gidhold)
test2 <- data.frame( track = stest$trackidhold, year = stest$yearhold, quad = stest$quadhold, Xhold  = stest$Xhold, group = stest$gidhold)

testdl <- merge( test1, test2, by = c('track', 'year', 'quad', 'group'))

plot(testdl$Xhold.x, testdl$Xhold.y) # scale is different but they match 
abline(0,1)

#
test1 <- data.frame( track = gtest$trackid, year = gtest$year, quad = gtest$quad, X  = gtest$X, group = gtest$gid)
test2 <- data.frame( track = stest$trackid, year = stest$year, quad = stest$quad, X  = stest$X, group = stest$gid)
testdl <- merge( test1, test2, by = c('track', 'year', 'quad', 'group'))
plot(testdl$X.x, testdl$X.y) # scale is different but they match 
abline(0,1) 

test1$Xu <- test1$X*gtest$Xscale + gtest$Xcenter
test2$Xu <- test2$X*stest$Xscale + stest$Xcenter

testdl <- merge( test1, test2, by = c('track' , 'year', 'quad', 'group'))
plot(testdl$Xu.x, testdl$Xu.y) # re-scaled
abline(0,1)

# 
test1 <- data.frame( track = gtest$trackid, year = gtest$year, quad = gtest$quad, X  = gtest$X, group = gtest$gid)
test2 <- data.frame( track = gtest$trackid3, year = gtest$year3, quad = gtest$quad3, X  = gtest$X3, group = gtest$gid3)
testdl <- merge( test1, test2, by = c('track', 'year', 'quad', 'group'))
plot(testdl$X.x, testdl$X.y)
abline(0,1) 
testdl[ testdl$X.x != testdl$X.y,  ]


# 
gdat <- readRDS('data/temp_data/ARTR_growth.RDS')
sdat <- readRDS('data/temp_data/ARTR_survival.RDS')
sdat <- sdat[ , c('quad', 'year', 'trackID', 'Group', 'logarea', 'survives', 'W.POSE')]
gdat <- gdat[ , c('quad', 'year', 'trackID', 'Group', 'logarea.t0', 'W.POSE')]

nrow(gdat)
nrow(sdat)

both <- merge(gdat, sdat, all = T, by = c('Group', 'quad', 'year', 'trackID'))
nrow( both[ is.na(both$logarea.t0), ] ) 
both[is.na(both$logarea.t0), ] 
plot(both$logarea, both$logarea.t0)
plot( both$W.POSE.x, both$W.POSE.y)
abline( 0, 1)
