df <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS') 

colnames(  df$ARTR$C  ) 
pairs( df$ARTR$C ) 
cor( unique( df$ARTR$C[ , 1:4]) ) 

colnames(  df$HECO$C  ) 
pairs( df$HECO$C ) 
cor( unique( df$HECO$C) ) 

colnames(  df$POSE$C  ) 
pairs( df$POSE$C ) 
cor( unique( df$POSE$C) ) 

colnames(  df$PSSP$C  ) 
pairs( df$PSSP$C ) 
cor( unique( df$PSSP$C[ , 1:3]) ) 

df <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS') 

colnames(  df$ARTR$C  ) 
pairs( df$ARTR$C ) 
cor( unique( df$ARTR$C[ , 1:2]) ) 

names(  df$HECO$Ccenter  ) 
# pairs( df$HECO$C ) 
# cor( unique( df$HECO$C) ) 

colnames(  df$POSE$C  ) 
pairs( df$POSE$C ) 
cor( unique( df$POSE$C [ , 1:3]) ) 

colnames(  df$PSSP$C  ) 
pairs( df$PSSP$C ) 
cor( unique( df$PSSP$C[ , 1:2]) ) 

df <- readRDS('data/temp_data/modified_recruitment_data_lists_for_stan.RDS') 

colnames(  df$ARTR$C  ) 
pairs( df$ARTR$C ) 
cor( unique( df$ARTR$C) ) 

names(  df$HECO$Ccenter  ) 
pairs( df$HECO$C ) 
cor( unique( df$HECO$C) ) 

colnames(  df$POSE$C  ) 
pairs( df$POSE$C ) 
cor( unique( df$POSE$C[ , 1:2] )  )

colnames(  df$PSSP$C  ) 
pairs( df$PSSP$C ) 
cor( unique( df$PSSP$C[ , 1:2]) ) 

