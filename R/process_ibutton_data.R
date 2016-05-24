rm(list = ls())

convert_time <- function(x) { strptime(x = x$Date.Time, format = '%m/%d/%y %I:%M:%S %p') } 

change_time <- function(x) { 
  x$date <- convert_time( x )
  return(x)
}

folders <- dir('data/iButton_data', recursive = FALSE , full.names = TRUE) 

for( i in 1:length(folders)){ 
  
  record_file <- dir( folders[i] , pattern = 'record', full.names = TRUE, recursive = TRUE)

  record <- read.csv(record_file)

  datafiles <- dir(folders[i], pattern = '[2|3].*21\\.csv', full.names = TRUE)

  d <- lapply( datafiles, read.csv, skip = 14)

  names(d) <- gsub(pattern = '.csv', replacement = '', basename(datafiles))

  d <- lapply(d, change_time)

  df <- do.call(rbind, d)

  df$id <- gsub( row.names(df), pattern = '\\.[0-9]+', replacement = '')

  df <- merge( df, record , by.x  = 'id', by.y  = 'ibutton')
  
  saveRDS(df, file = file.path(folders[i], 'processed.RDS') )
}



