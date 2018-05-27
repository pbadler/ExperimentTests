# 
library(stringr)
library(dplyr)
library(tidyr)

rm(list = ls())

cor_files <- dir('output', 'correlations.csv', full.names = T)
cor_files <- cor_files [ -grep('all', cor_files)] # keep only the selected ones 

out <- out2 <-  list( NA )
i = 1
for(i in 1:length(cor_files)){
  
  spp <- str_extract( pattern =  '[A-Z]{4}', cor_files[i])
  vr  <- str_extract( pattern = '(growth)|(recruitment)|(survival)', cor_files[i])
  
  temp <- read.csv(cor_files[i])
  
  out2[[i]] <- temp
  out2[[i]]$species <- spp
  out2[[i]]$vital_rate <- vr
  names(out2)[i] <- paste(vr, spp)
  
  
  VWC_parameter <- as.character( temp$var [ temp$vartype == 'VWC'] )
  
  ############################################
  # hard code which species recieve size by climate interaction 
  # based on whether those species have size by treatment interaction 

  if (vr == 'growth' & spp == 'PSSP') {  # use the size effect interaction 
    ifx <- paste( VWC_parameter, 'logarea.t0', sep = 'x')
    
    VWC_parameter <- c(VWC_parameter, ifx )
  }

  if( vr == 'survival' & spp == 'POSE' ) { # use the size effect interaction 
    ifx <- paste(VWC_parameter, 'logarea.t0', sep = 'x')
    VWC_parameter <- c(VWC_parameter, ifx)
  }  
  #############################################
  
  use_vars <- VWC_parameter
  
  out[[i]] <- data.frame( species = spp, vital_rate = vr, covars = paste( use_vars, collapse = ','))
  
  
}

out <- do.call(rbind, out)

write.csv(out, 'output/selected_climate_covariates.csv')

#

recOuts <- out2[grep('recruitment', out2)]
gsOuts <- out2[-grep('recruitment', out2)]

out_cors <- lapply( gsOuts, function(x) x[, c('Intercept', 'Int_pval', 'size', 'size_pval', 'var', 'species', 'vital_rate')])

out_cors <- do.call(rbind, out_cors)

rec_outs <- do.call(rbind, recOuts)


rec_outs$size <- NA
rec_outs$size_pval <- NA
rec_outs <- rec_outs[, c('Intercept', 'Int_pval', 'size', 'size_pval', 'var', 'species', 'vital_rate')]

all_cors <- rbind(out_cors, rec_outs)

all_cors <- all_cors %>% dplyr::select(vital_rate, species, var, Intercept, Int_pval, size, size_pval ) %>% arrange( vital_rate, species, desc(abs(Intercept)) )

################################################
# hard code which species recieve size by climate interaction 
# based on whether those species have size by treatment interaction 

size_cors <- 
  all_cors %>% 
  filter( ((species == 'PSSP' & vital_rate == 'growth') | ( species == 'POSE'  & vital_rate == 'survival')) )

all_cors <- 
  all_cors %>% 
  filter( !((species == 'PSSP' & vital_rate == 'growth') | ( species == 'POSE'  & vital_rate == 'survival')))

################################################

all_cors$size <- NA
all_cors$size_pval <- NA

all_cors <- rbind( all_cors, size_cors)

all_cors <- all_cors %>% rename(`climate variable` = var , `Int. cor.` = Intercept, `p val.`  = Int_pval, `Size cor.` = size , `Size p. val.`  = size_pval )

all_cors <- all_cors %>% arrange( vital_rate, species, desc(abs(`Int. cor.`)))

all_cors$`climate variable` <- str_replace( all_cors$`climate variable` , 'l', 'lag')
all_cors$`climate variable` <- str_replace(all_cors$`climate variable`, 'VWC.', '')


library(xtable)


xtcor <- xtable(all_cors, 
                caption = 'Selected climate variables for each vital rate model for each species. Correlations and p-values between the choosen variables and the intercept of the no climate model are shown. For ARTR growth and POSE growth and survival, the correlations between the year effects on size and the soil moisture variables are also given. "f" = fall, "su" = summer, "sp" = spring. ARTR = \\textit{A. tripartita}, HECO = \\textit{H. comata}, POSE = \\textit{P. secunda}, PSSP = \\textit{P. spicata}.', 
                label = 'table:strongCor')

print(xtcor, 'manuscript/strong_correlations.tex', type = 'latex', include.rownames = F, caption.placement ="top")

# check for highly correlated covariates 
# clim <- readRDS('data/temp_data/all_clim_covs.RDS')
# clim <- clim[complete.cases(clim), ]
# correlations <- list()
# i = 3
# for ( i in 1:nrow( out )  ) { 
#   covars     <- strsplit( as.character( out$covars[i] ) , ',')[[1]]
#   covars <- covars[ -grep('x', covars)]
#   correlations[[i]] <- cor(clim[ , covars, drop = F] )
# }
# 
# all_cors <- lapply(cor_files,  read.csv)
# 
# names(all_cors) <- basename( cor_files )
# 
# recruitment_cors <- do.call( rbind, all_cors[ grep( 'recruitment', names(all_cors)) ] )
# 
# labels <- do.call( rbind, str_split(row.names(recruitment_cors), '_'))[, 1:2]
# recruitment_cors$species <- labels[, 1]
# recruitment_cors$vital_rate <- labels[, 2]
# 
# gs_cors <- do.call( rbind, all_cors[ -grep( 'recruitment', names(all_cors)) ] )
# 
# labels <- do.call( rbind, str_split(row.names(gs_cors), '_'))[, 1:2]
# gs_cors$species <- labels[, 1]
# gs_cors$vital_rate <- labels[, 2]
# 
# recruitment_cors
# gs_cors
# 
# write.csv(file = 'output/year_effects_correlations_recruitment.csv', recruitment_cors, row.names = F)
# write.csv(file = 'output/year_effects_correlations_growth_survival.csv', gs_cors, row.names = F)
# 


