#### plot daily max temps 
rm(list = ls())

library(ggplot2)
source( 'retrieveMesoWestData.R')

lastDay = Tmax[which.max(Tmax$doy), ]
ylab = 'Max Daily Temp in Celsius'
xlab = 'Day of year'

plotTitle = paste( 'Daily Tmax at US Sheep Experiment Station\n 
              from 1925-01-01 to', strftime(lastDay$date) )

theme_set( theme_bw() )

p1 = ggplot(dailyRecords, aes( x = doy, y= tmax, group = year)) +   
  geom_line(color = 'darkgray') + geom_line(data = subset( dailyRecords, year == 2015), color = 'purple', lwd = 1.5) + 
  geom_line( data = TmaxStats, aes( x = doy, y = TmaxMN, group = 1), color = 'red') + 
  geom_line( data = TmaxStats, aes( x = doy, y = TmaxMN + 2*TmaxSD, group = 1), color = 'red') + 
  geom_line( data = TmaxStats, aes( x = doy, y = TmaxMN - 2*TmaxSD, group = 1), color = 'red') + 
  ylab( ylab ) + xlab( xlab ) + ggtitle(plotTitle) 

#p1 = p1 + geom_line( data = subset(dailyRecords, V1 == '2013'), aes( x = doy, y = V4, group = 1), color = 'darkgreen', lwd = 1.5)

###### Draw on MesoWest recent data 
p2 = p1 + geom_line( data = subset(Tmax, year == 2015), aes( x = doy, y = temp_C), color = 'purple', lwd = 1.5) + 
  geom_point( data = lastDay, aes(x = doy, y = temp_C), color = 'purple') +
  geom_text( data = lastDay, aes(x = doy, y = temp_C, label = year), color = 'purple', hjust = -0.1) + 
  xlim( c(0, lastDay$doy + 31 )) + geom_hline( yintercept = 0, alpha = 0.3)

p2

###### Find number of days with Tmax > average Tmax 
TmaxComparison = merge( Tmax, TmaxStats, by = 'doy', suffixes=c('.avg', '.2015'))
head( TmaxComparison)

head( Tmax ) 
TmaxStats

older = dailyRecords[ which( !dailyRecords$doy %in% Tmax$doy & dailyRecords$year == 2015), c('tmax', 'doy')]
newer = Tmax[ , c('temp_C', 'doy')] 
names(older) <- c('Tmax2015', 'doy')
names(newer) <- c('Tmax2015', 'doy')

curYear = rbind(older, newer)
TmaxComparison = merge( curYear, TmaxStats, by = 'doy')

TmaxComparison$compare <- NA 
sum( TmaxComparison$Tmax2015 > TmaxComparison$TmaxMN) 
sum(TmaxComparison$Tmax2015 < TmaxComparison$TmaxMN)

length(which(TmaxComparison$Tmax2015 > TmaxComparison$TmaxHi))
length(which(TmaxComparison$Tmax2015 < TmaxComparison$TmaxMN ))

TAnom2015 =  TmaxComparison$Tmax2015 - TmaxComparison$TmaxMN  #### Temperature Anomaly 
Tvar = TmaxComparison$TmaxHi - TmaxComparison$TmaxLo
plot(Tvar, type = 'l') 

mean( TAnom2015 ) 
plot(TAnom2015, type = 'l')
abline(h = 0 )
hist(TAnom2015)
abline( v= mean(TAnom2015))

p2 + geom_smooth(data = curYear[ curYear$doy < 300, ], aes( x = doy, y = Tmax2015), color = 'blue', linetype = 2)
