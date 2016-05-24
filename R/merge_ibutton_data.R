rm(list = ls()) 

q_info <- read.csv('data/quad_info.csv') 

f <- dir('data/iButton_data', pattern = '.*\\.RDS', recursive = TRUE, full.names = TRUE)

df <- lapply( f, readRDS )

df <- do.call(rbind, df )

head(q_info)

q_info$plot <- gsub( q_info$QuadName, pattern = 'X', replacement = '')

df <- merge( df, q_info, by = 'plot')

ggplot( df , aes( x = date, y = Value, color = Treatment , group = QuadName )) + 
  geom_point() + 
  geom_line() + facet_wrap( ~ PrecipGroup, ncol = 1)


