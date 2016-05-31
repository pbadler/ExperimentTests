library('ggplot2')
library('SPEI')
library('xtable') 


setwd('~/Desktop/sheep_station_historical_climate/')


daily = read.csv('Dubois_Climate_History.csv')

daily$Date = as.POSIXct(strptime( daily$DATE, "%Y%m%d"))

daily$DOY = strftime(daily$Date, format = "%j")
daily$year = strftime(daily$Date, format = "%Y")
head(daily)

yearlyTotals = aggregate(PRCP ~ year, subset(daily, PRCP != -9999), FUN = 'sum')

recent = subset(daily, year > 2000 & year < 2014 )
recent$WaterYear = NA
recent$WaterYear = c(rep(2001, 273), recent$year[ which( recent$year >  2001 ) ], rep(2014, 92) ) 
recent <- recent[ recent$WaterYear > 2001 & recent$WaterYear < 2014 , ]

recent$DOWY = 1 

recent$DOWY = as.numeric(unlist(aggregate( DOWY ~ WaterYear, recent, 'cumsum')))[ -c(1:12)]

recent[recent$PRCP < 0 , "PRCP" ] <- 0
csum = aggregate(PRCP ~ WaterYear, recent, FUN = 'cumsum')
recent$cumsum = unlist(csum$PRCP)

totals = aggregate(PRCP ~ WaterYear, recent, FUN = 'sum')

ggplot(recent, aes( x = as.numeric(DOWY), y = cumsum, color = WaterYear)) + geom_line()
totals$DM50 = NA
head(recent)
for ( i in 1:nrow(totals)){ 
  year = totals$WaterYear[i]
  cutoff = totals$PRCP[i]*0.5
  totals$DM50[i] = max( recent$DOWY[ which(recent$WaterYear == year & recent$cumsum <= cutoff ) ] )
  totals$date50[i] = recent$DATE[ which(recent$WaterYear == year & recent$DOWY == totals$DM50[i]) ]
}


totals

ggplot(totals, aes(x = PRCP, y = DM50)) + geom_point() + geom_text( aes( label = WaterYear), adj = 1)

totals [ - c(1, 14), ]
