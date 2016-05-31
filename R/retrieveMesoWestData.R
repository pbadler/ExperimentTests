
require(RCurl)
require(RHTMLForms)
require(XML)

source('getDailyMaxTemps.R')

loginurl = "http://mesowest.utah.edu/"
pword = 'ring218'
user = 'arklein'

baseURL = paste( loginurl, 'cgi-bin/droman/meso_download_mesowest_ndb.cgi?product=', sep = '' )
##### example options: &stn=DUB&unit=0&time=LOCAL&daycalendar=1&day1=18&month1=03&year1=2015&yearcal=2015&hour1=10&monthcal=03&hours=90&output=csv&order=1&TMPF=TMPF

overlap = 20  #### num days of overlap between most recent mesowest and climate explorer data

#### Look-up parameters 
yesterday = Sys.Date() - 1
day = strftime( yesterday, format = '%d') 
month = strftime(yesterday, format = '%m')
year = strftime( yesterday, format = '%Y')
ndays = as.numeric( yesterday - as.Date(lastDay) ) + overlap  ##### number of days to retrieve

stationName = 'DUB'
timeOpt = 'LOCAL'
daycalendar = 1
hour1 = 00
output = 'csv'
orderOpt = 1
vars = 'TMPF'
unit = 1

loginForm = getHTMLFormDescription(loginurl) #### read website and extract all forms information
loginFun = createFunction(loginForm[[2]]) ##### create login function to pass information to the form 

dataURL = paste( baseURL, '&stn=', stationName, '&unit=' , unit, '&time=' , timeOpt, '&daycalendar=', daycalendar, '&day1=', day, '&month1=', month, '&year1=', year,
                 '&yearcal=', year, '&hour1=', hour1, '&monthcal=', month, '&hours=', ndays, '&output=', output, '&order=', orderOpt, '&', vars, '=', vars, sep = '') 

#### sets up curl 
agent="Mozilla/5.0" 
curl = getCurlHandle()
curlSetOpt(cookiefile = "", useragent = agent, verbose = TRUE, followlocation = TRUE, curl=curl)

loginFun(password= pword, user= user, .curl= curl)

##### parse data 
data = getURLContent( dataURL, curl = curl)

data.parsed = htmlTreeParse(data, useInternal = TRUE)
root = xmlRoot(data.parsed)
text = xmlValue(root)

pattern = '[A-Z]+[\n]{2}([0-9].*)'
match = regexec( pattern= pattern, text= text)
datTable = regmatches( x= text, m= match )[[1]][2]
datTable = textConnection( datTable )

datDF = read.csv(datTable, header = FALSE)
names(datDF) <- c('month', 'day', 'Y', 'hour', 'min', 'tz', 'temp_C')

##### Clean data 
datDF$date = strptime( x= paste( datDF$Y, datDF$month, datDF$day, sep = '-'), format = '%Y-%m-%d' )
datDF$doy = as.numeric(strftime( datDF$date, format = '%j'))

Tmax = aggregate( temp_C ~  Y + month + day, datDF, FUN= 'max')
Tmax$date = strptime( x= paste( Tmax$Y, Tmax$month, Tmax$day, sep = '-'), format = '%Y-%m-%d' )
Tmax$doy = as.numeric(strftime( Tmax$date, format = '%j'))
Tmax$year = as.numeric(strftime(Tmax$date, format = '%Y'))

#saveRDS( object= Tmax, file='recentMaxTempsMesoWest.RDS')

