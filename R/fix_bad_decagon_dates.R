# make file specific corrections 

# trouble files: 

# plot EM20086 period 3: 
#   12 AM and 12 PM are out of order in the sheet so the readings are out of order 
# 
# plot EM6843 period 3: 
#   dates jump up to 2057 in the data sheet. then jump back to 2013. 
#   could correct by assigning the dates > 2013 consecutive dates starting with the last 2013 date before the jump. 
#
# plot 11_12_C period 4:  
#   skips from 11-09-2013 to 11-12-2013
#   jumps backwards 6 hours on 03-12-2014
#
# plot 12 period 4: 
#   skips 6 hours on 09-03-2013 
#   skips 12 hours on 09-05-2013 
#   skips 10 hours on 09-04-2013
#   skips from 11-09-2013 to 11-12-2013
#   skips 4 hours on 03-12-2014
#   skips 6 hours on 11-12-2013
# 
# plot 15_16_C period 4: 
#   skips back 2 hours and makes two separate readings at 4pm on 11-09-2013
#
# plot 7_8_C period 4: 
#   skips up 4 hours on 11-09-2013 
#   jumps back 6 hours on 03-12-2014
#
# plot 8 period 4: 
#   skips up 4 hours on 11-09-2013 
#   jumps back 6 hours on 03-12-2014
#
# plot 15 , period 5: 15_15Sep14-1802.txt
#   readings are out of order with respect to date, but dates look correct 
# 
# plot EM20068 period 6: 
#   jumps from 03-28-2015 to 08-03-2015 
#   fix by adjusting dates > 08-02-2015
#
# plot EM20085 period 6: 
#   jumps back from 10-27-2014 to 01-02-2000
#   fix by adjusting dates < 01-01-2011
#
# plot EM20086 period 6: 
#   skips hours in 03-03-2015
#   skips hours on 03-03-2015
#   skips hours on 03-04-2015
#   skips hours on 03-18-2015
#   skips hours on 03-19-2015
#   skips hours on 03-20-2015
#
# plot EL5739 period 7:
#   skips hours in July and Aug 2015
#
# plot EL5742 period 7: 
#   skips hours on 09-18-2015
#
# plot EL5743 period 7:
#   skips hours in July and Aug 2015
#
# plot EM20068 period 7:
#   jumps from 03-28-2015 to 08-03-2015  (corrected for period 6)
#   jumps from 09-05-2015 to 04-30-2015 
#     fix by adjusting 
#   jumps from 06-27-2015 to 02-01-2049 
#     fix by adjusting 
#   jumps from 02-04-2049 to 01-02-2000
#     fix by adjusting 


convert_time <- function(x) { 
  
  return( strptime(x = x$Time, format = '%m/%d/%y %I:%M %p', tz = 'MST') ) 
  
}


# ---------------------------------------------------------------------------------------------------------------
# correction to EM6843 period 3 
# file: EM6843_2013_11_09-processed.txt'
# ---------------------------------------------------------------------------------------------------------------

EM6843 <- read.table( 'data/soil_moist_data/2013_Fall/EM6843_2013_11_09-processed.txt', sep = '\t', skip = 1 )

EM6843$Time <- EM6843$V1
date <- convert_time ( EM6843 )

tdiff <- diff.POSIXt(date)

table ( as.numeric( tdiff, units = 'secs') ) # use the larger difference here as the correction 

correction <- 1398427200 - 2*60*60 # correction, converted to seconds, offset by two hours  

df <- readRDS(file = 'data/temp_data/decagon_data.RDS')

df <- df %>% 
  mutate( new_date = ifelse( id == 'EM6843' & 
                               period == 3 & 
                               date > as.POSIXct( '2014-01-01', '%Y-%m-%d', tz = 'MST'), 
                             as.POSIXct(date - correction, origin = '1970-01-01 00:00:00', tz = 'MST'), # if yes 
                             new_date ))                                                                # if no  

test <- df %>% filter( id == 'EM6843', period == 3 )
test <- test %>% arrange ( reading )
plot ( data = subset( test , depth == 'air' ), value~ new_date )  # check 
plot ( data = subset( test, depth == 'air'), value  ~ reading ) # check 

saveRDS(df, file = 'data/temp_data/decagon_data.RDS') # replace file 

# ---------------------------------------------------------------------------------------------------------------
# correction to plot EM20068  period 6 
# file: /home/andy/Documents/precip_experiment/data/soil_moist_data/2015_Spring/EM20068_2015-04-30-0957.txt
# ---------------------------------------------------------------------------------------------------------------

EM20068 <- read.table( 'data/soil_moist_data/2015_Spring/EM20068_2015-04-30-0957.txt', sep = '\t', skip = 1 )

EM20068$Time <- EM20068$V1
date <- convert_time ( EM20068 )

tdiff <- diff.POSIXt(date)

table ( as.numeric( tdiff, units = 'secs') ) # use the larger difference here as the correction 

correction <- 11052000 - 2*60*60 # correction, converted to seconds, offset by two hours  

df <- readRDS(file = 'data/temp_data/decagon_data.RDS')

df <- df %>% 
  mutate( new_date = ifelse( id == 'EM20068' & 
                               period == 6 & 
                               date > as.POSIXct( '2015-08-01', '%Y-%m-%d', tz = 'MST'), 
                             as.POSIXct(date - correction, origin = '1970-01-01 00:00:00', tz = 'MST'), 
                             new_date )) 

test <- df %>% filter( id == 'EM20068', period == 6 )
test <- test %>% arrange ( reading )
plot ( data = subset( test , depth == 'air' ), value~ new_date )  # check 
plot ( data = subset( test, depth == 'air'), value  ~ reading ) # check 
plot ( data = subset( test, depth = 'air'), value ~ date )

saveRDS(df, file = 'data/temp_data/decagon_data.RDS') # replace file 

# ---------------------------------------------------------------------------------------------------------------
# correction to plot EM20085  period 6 
# file: 2015_Spring/EM20085_2015-04-30-1008.txt
# ---------------------------------------------------------------------------------------------------------------

EM20085 <- read.table( 'data/soil_moist_data/2015_Spring/EM20085_2015-04-30-1008.txt', sep = '\t', skip = 1 )

EM20085$Time <- EM20085$V1
date <- convert_time ( EM20085 )

tdiff <- diff.POSIXt(date)

table ( as.numeric( tdiff, units = 'secs') ) # use the larger difference here as the correction 

correction <- -467611200 - 2*60*60 # correction, converted to seconds, offset by two hours  

df <- readRDS(file = 'data/temp_data/decagon_data.RDS')

df <- df %>% 
  mutate( new_date = ifelse( id == 'EM20085' & 
                               period == 6 & 
                               date < as.POSIXct( '2011-01-01', '%Y-%m-%d', tz = 'MST'), 
                             as.POSIXct(date - correction, origin = '1970-01-01 00:00:00', tz = 'MST'), 
                             new_date )) 

test <- df %>% filter( id == 'EM20085', period == 6 )
test <- test %>% arrange ( reading )
plot ( data = subset( test , depth == 'air' ), value~ new_date )  # check 
plot ( data = subset( test, depth == 'air'), value  ~ reading ) # check 
plot ( data = subset( test, depth = 'air'), value ~ date )

saveRDS(df, file = 'data/temp_data/decagon_data.RDS') # replace file 

# ---------------------------------------------------------------------------------------------------------------
# correction to plot EM20068 period 7 
# file: 2015_Fall/EM20068
# ---------------------------------------------------------------------------------------------------------------

EM20068 <- read.table( 'data/soil_moist_data/2015_Fall/EM20068 4Nov15-1138.txt', sep = '\t', skip = 1 )

EM20068$Time <- EM20068$V1
date <- convert_time ( EM20068 )

tdiff <- diff.POSIXt(date)

jumps <- data.frame(tdiff = c(NA, as.numeric(tdiff, units = 'secs')), date, reading = 1:nrow(EM20068) ) %>% distinct ( tdiff ) 

jumps # view time jumps 

jumps$correction <- jumps$tdiff - 2*60*60 # corrections are jumps converted to seconds and offset by two hours 

df <- readRDS(file = 'data/temp_data/decagon_data.RDS')

test <- df %>% filter( id == 'EM20068', period == 7 )

df <- df %>% 
  mutate( new_date = ifelse( id == 'EM20068' & 
                               period == 7 & 
                               reading >= jumps$reading[3],
                             as.POSIXct(date - jumps$correction[3], origin = '1970-01-01 00:00:00', tz = 'MST'), new_date ))

df <- df %>% 
  mutate( new_date = ifelse( id == 'EM20068' &
                               period == 7 &
                               reading > jumps$reading[4], 
                             as.POSIXct(date - jumps$correction[4], origin = '1970-01-01 00:00:00', tz = 'MST'), new_date ))
df <- df %>% 
  mutate( new_date = ifelse( id == 'EM20068' &
                               period == 7 &
                               reading > jumps$reading[5], 
                             as.POSIXct(date - jumps$correction[5], origin = '1970-01-01 00:00:00', tz = 'MST'), new_date ))

df <- df %>% 
  mutate( new_date = ifelse( id == 'EM20068' &
                               period == 7 &
                               reading > jumps$reading[6], 
                             as.POSIXct(date - jumps$correction[6], origin = '1970-01-01 00:00:00', tz = 'MST'), new_date))

test <- df %>% filter( id == 'EM20068', period == 7 )
test <- test %>% arrange ( reading )
plot ( data = subset( test , depth == 'air' ), value~ new_date )  # check 
plot ( data = subset( test, depth == 'air'), value  ~ reading ) # check 
plot ( data = subset( test, depth = 'air'), value ~ date )

saveRDS(df, file = 'data/temp_data/decagon_data.RDS') # replace file 

