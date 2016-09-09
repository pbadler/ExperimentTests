rm(list = ls())

q_info <- read.csv('data/quad_info.csv') 

folders <- dir('data/iButton_data', recursive = FALSE , full.names = TRUE) 

data_list <- list(NA)

for( i in 1:length(folders)){  # process each folder 
  
  record_file <- dir( folders[i] , pattern = 'record', full.names = TRUE, recursive = TRUE) 

  record <- read.csv(record_file)

  datafiles <- dir(folders[i], pattern = '[2|3].*21\\.csv', full.names = TRUE) # list all data files in the folder 

  d <- lapply( datafiles, read.csv, skip = 14)  

  names(d) <- gsub(pattern = '.csv', replacement = '', basename(datafiles))

  df <- do.call(rbind, d)

  df$date <- strptime(df$Date.Time, format = '%m/%d/%y %I:%M:%S %p' )
  
  df$date <- as.POSIXct(df$date)
  
  df$id <- gsub( row.names(df), pattern = '\\.[0-9]+', replacement = '')

  df <- merge( df, record , by.x  = 'id', by.y  = 'ibutton')
  
  data_list[[i]] <- df 
}


df <- do.call( rbind, data_list )  # bind the data lists from each folder 

q_info$plot <- gsub( q_info$QuadName, pattern = 'X', replacement = '')

df <- merge( df, q_info, by = 'plot') 

saveRDS(df, 'data/temp_data/ibutton_data.RDS')

