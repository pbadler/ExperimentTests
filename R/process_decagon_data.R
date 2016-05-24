rm(list = ls())

library(stringr)
library(dplyr)
library(tidyr)

get_col_names <- function( x ) { 
  
  port <- str_extract_all(x[1, ], pattern = '(Time)|(Port\\s[0-9])') 
  probe_type <- str_extract( x[1, ] , pattern = '(ECT)|(EC\\-5)|(5TM)')
  measure <- str_extract( x [ 1, ] , pattern = '(VWC$)|((C$|C\\sTemp$))')
  
  new_names <- paste( port, measure , sep = '_')
  str_replace_all(string = new_names, pattern = c('_NA'), replacement = c('')) 
  
}


rename_cols <- function(x ) { 
  
  names(x) <- get_col_names( x ) 
  
  return(   x[-1, ] ) 
}

assign_NAs <- function( x ) { 
  
  x [ x == '* * * '] <- NA 
  
  return( x )
}

convert_time <- function(x) { 
  
  strptime(x = x$Time, format = '%m/%d/%y %I:%M %p') 
}


change_time <- function(x) { 
  
  x$date <- convert_time( x )
  
  x$date <- as.POSIXct(x$date )
  
  return(x)
}

gather_ports <- function ( test ) { 
  test %>% 
    gather( key = port, value = value ,  starts_with("Port") ) %>% 
    separate(col = port, into = c('port', 'measure') , sep = '_') 
} 

q_info <- read.csv('data/quad_info.csv') 

folders <- dir('data/soil_moist_data', pattern = '20[0-9]{2}_((Fall$)|(Spring$))', full.names = TRUE)

data_list <- list(NA)

for (i in 1:length(folders)) {
  i = 1
  record_file <- dir(folders[i] , pattern = 'logger_info.csv', full.names = TRUE) 

  record <- read.csv(record_file)
  
  f <- dir(folders[i], pattern = '^E[ML][0-9]+.*txt$', full.names = TRUE)
  
  f2 <- dir(folders[i], pattern = '^[0-9]+(_[0-9]+_C)?.*txt$', full.names = TRUE)
  
  f <- c(f, f2)
  
  d <- lapply(f, read.table, sep = '\t', colClasses = 'character')  
  
  names(d) <- str_extract(basename(f), pattern = '(^E[ML][0-9]+)|(^[0-9]+(_[0-9]+_C)?)')

  d <- lapply(d, rename_cols) 

  d <- lapply(d, assign_NAs)

  d <- lapply(d, change_time) 

  d <- lapply( d, gather_ports ) 

  df <- do.call(rbind, d)

  df$id <- gsub( pattern = '\\.[0-9]+$', replacement = '', x = row.names(df))

  df <- merge(df, record, by.x = 'id', by.y = 'logger' )
  
  df$value <- as.numeric(df$value)
    
  data_list[[i]] <- df 
  
} 


df <- do.call( rbind, data_list )  # bind the data lists from each folder 

q_info$plot <- gsub( q_info$QuadName, pattern = 'X', replacement = '')

df <- merge( df, q_info, by = 'plot') 

saveRDS(df, 'data/temp_data/decagon_data.RDS')