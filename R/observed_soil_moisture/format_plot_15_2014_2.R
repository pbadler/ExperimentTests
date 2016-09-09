# special formatting script for plot 15 

df <- read.table(file = 'data/soil_moist_data/2014_2/15_15Sep14-1802.txt', sep = '\t', header = F)

cn <- df[ 1, ]

df <- df [ -1 , ]

df$date <- strptime( df$V1, '%m/%d/%y %I:%M %p', tz = 'MST' ) 

df <- df[ order(df$date), ]

df <- df [ , -7 ]

df <- rbind( cn, df )

write.table( df, 'data/soil_moist_data/2014_2/15_reorderd.txt', sep = '\t', row.names = F, col.names = F, quote = F)
