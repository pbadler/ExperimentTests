library('ggplot2')
library('SPEI')
library('xtable') 


setwd('~/Desktop/sheep_station_historical_climate/')

Lat = 44.15 #### station latitude

daily = read.csv('Dubois_Climate_History.csv')
clim = read.csv('Dubois_Climate_monthly.csv')
clim$Date = as.POSIXct(strptime( clim$DATE, "%Y%m%d"))
head(clim, 13)
tail(clim, 32)
clim[ clim == -9999] = NA

clim$year = as.numeric(strftime(clim$Date, '%Y'))
clim$month = as.numeric(strftime(clim$Date, '%m'))
tail(clim, 10)

clim$wYear = c(clim$year[ -c(1:3) ], c(2014, 2014, 2014))
clim$seasonal = FALSE
clim$seasonal[clim$month %in% c(1,2,3,4,10,11,12) ] <- TRUE
tail(clim, 10)

PETmat = hargreaves(Tmin=clim$MMNT, Tmax= clim$MMXT, Pre= clim$pptmm, lat = Lat, na.rm = TRUE)

waterClim = data.frame(year = clim$year, wYear = clim$wYear, month = clim$month, data.frame(PETmat), PPTmm = clim$TPCP/10, wetSeason = clim$seasonal)

head(waterClim)
seasonalPET = aggregate(ET0_har ~ wYear, data = subset(waterClim, wetSeason == TRUE), FUN = 'sum') 
seasonalPPT = aggregate(PPTmm ~ wYear, data = subset(waterClim, wetSeason == TRUE), FUN = 'sum') 
seasonalDEF = data.frame(DEF = seasonalPET$ET0_har - seasonalPPT$PPTmm, wYear = seasonalPET$wYear)

dryRanked = seasonalPPT[ order(seasonalPPT$PPTmm), ]
dryest = dryRanked[1:10, ]
wetest = tail(dryRanked, 10)
dryRanked$recent <- FALSE
dryRanked$recent [ dryRanked$wYear > 2006 ] <- TRUE
PETRanked = seasonalPET[ order(seasonalPET$ET0_har, decreasing= TRUE), ]
PEThigh = PETRanked[1:10, ]
DEFRanked = seasonalDEF[ order(seasonalDEF$DEF, decreasing= TRUE), ]
DEFhigh = DEFRanked[1:10, ]


seasonalPPT
labeldf = data.frame(wYears = c(2012, 2013, 2014), vals = tail(seasonalPPT$PPTmm, 3))
labeldf

theme_set(theme_classic(base_size= 10)) 

histPlot = ggplot(aes(), data = seasonalPPT ) + geom_histogram( aes(x = PPTmm), bin = 20) + 
  xlab('Oct - April precip (mm)') + ylab('Number of years') + 
  geom_text( data = labeldf, aes(label = wYears, x = vals, y = -1), col = 'black', size = 3.5) +
  geom_vline(x = labeldf$vals, col = "red", lty = 2) 


ENSO = data.frame( type = c(rep('ElNino', 6), rep('LaNina', 5)), year = c(1958, 1966, 1973, 1983, 1988, 1998, 1974, 1976, 1989, 2000, 2011))

#par(mfrow = c(2, 1))
#plot(seasonalPPT$wYear, seasonalPPT$PPTmm, type = 'l', xlab = '', ylab = "Total Precip. (mm)")
#abline(v = ENSO$year[ ENSO$type == 'ElNino'] , col = 'red')
#abline(v = ENSO$year[ ENSO$type == 'LaNina'], col = 'blue')
#text(label = dryest$wYear, x = dryest$wYear, y = dryest$PPTmm)
#plot(seasonalPET$wYear, seasonalPET$ET0_har, type = 'l', xlab = '', ylab = "Potential E.T. (mm)")
#text(label = PEThigh$wYear, x = PEThigh$wYear, y = PEThigh$ET0_har)


write.csv(dryRanked, "USSES_dry_years.csv", row.names = FALSE)

