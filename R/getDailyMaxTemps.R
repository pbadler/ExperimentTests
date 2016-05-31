#### 
#### Climate explorer data for Dubois, ID 
require(RCurl)
require(XML)

agent="Mozilla/5.0" 
curl = getCurlHandle()
curlSetOpt(cookiefile = "", useragent = agent, verbose = TRUE, followlocation = TRUE, curl=curl)

dataURL = 'http://climexp.knmi.nl/data/xgdcnUSC00102707.dat'
data = getURLContent(dataURL, curl= curl)
datTable = textConnection( data )

dailyRecords = read.table(datTable, header = FALSE)
names(dailyRecords) <- c('year', 'month', 'day', 'tmax')

dailyRecords$dateChar = paste( dailyRecords$year, dailyRecords$month, dailyRecords$day, sep = '-')
dailyRecords$date = strptime( dailyRecords$dateChar, format = '%Y-%m-%d')
dailyRecords$doy = as.numeric( strftime( dailyRecords$date, '%j'))

TmaxSD = aggregate(tmax ~ doy, dailyRecords, FUN = 'sd')
TmaxMN = aggregate(tmax ~ doy, dailyRecords, FUN = 'mean' )
TmaxHi = aggregate(tmax ~ doy, dailyRecords, FUN = 'max')
TmaxLo = aggregate(tmax ~ doy, dailyRecords, FUN = 'min')

TmaxStats = cbind( TmaxMN, TmaxHi[,2], TmaxLo[,2], TmaxSD[,2])

names( TmaxStats ) <- c('doy', 'TmaxMN', 'TmaxHi', 'TmaxLo', 'TmaxSD')

dailyRecords = dailyRecords [ -nrow(dailyRecords), ] #### skip last day

lastDay = max(dailyRecords$date)

