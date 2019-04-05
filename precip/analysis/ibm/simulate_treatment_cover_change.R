rm(list =ls())
library(tidyverse)

# functions ---------------------------------------------------------------------------- 

get_recruit_area <- function(spp) { 
  
  dat <- read.csv(paste0('~/Dropbox/driversdata/data/idaho/speciesData/', as.character(spp), '/recSize.csv'))
  
  dat$area

}

#
load('analysis/figure_scripts/my_plotting_theme.Rdata')

quads <- read_csv('data/quad_info.csv')

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

fns <- dir('output/stan_fits', '.RDS', full.names = T)

fn_df <- 
  data.frame(fns) %>% 
  mutate( bname = basename(as.character(fns))) %>% 
  separate( bname, c('spp', 'vr', 'window', 'type'), sep = '_') %>% 
  mutate( type = str_remove(type, pattern = '\\.RDS$' )) 
  
fn_df <- 
  fn_df %>% 
  group_by( spp, vr , type ) %>% 
  arrange( window) %>% 
  mutate( rank = row_number() ) %>% 
  mutate( climate = rank == 1) %>% 
  select(-rank ) 

model_list <- 
  rbind( 
  fn_df %>% 
    filter( climate ) %>% 
    unite( type, vr, type, sep = '_') %>% 
    select(-window) %>% 
    spread( type, fns )
  , 
  fn_df %>% 
    select( - climate ) %>% 
    filter( window == 'none') %>% 
    unite( type , vr, type , sep = '_' ) %>% 
    select( -window) %>% 
    spread( type, fns) %>% 
    mutate( climate = F)
)


model_list <- model_list %>% ungroup()

# loop species and climate / non-climate IBMs 

i = 1
spp <- species_list[i]

gd <- readRDS( as.character(model_list$growth_data[1]))
rd <- readRDS( as.character(model_list$recruitment_data[1]))
sd <- readRDS( as.character(model_list$survival_data[1]))

gfit <- readRDS( as.character(model_list$growth_model[1]))
sfit <- readRDS( as.character(model_list$survival_model[1]))
rfit <- readRDS( as.character(model_list$recruitment_model[1]))

# generate predicted area per quad per year -------------- # 
S <- binomial(link='logit')$linkinv(rstan::extract( sfit, 'IBM_mu')$IBM_mu)
G <- rstan::extract( gfit, 'IBM_Y_hat')$IBM_Y_hat

# Test that the predictions can be rescaled correctly to cm^2 #  
Y_attrib <- gd$IBM_Y_attrib
Y_center <- Y_attrib$`scaled:center`
Y_scale <- Y_attrib$`scaled:scale`

gdat <- readRDS('data/temp_data/ARTR_growth_survival_dataframe.RDS')

Y1 <- (gd$IBM_Y*Y_scale + Y_center)
Y2 <- gdat$logarea.t1
all.equal(Y1, Y2)
#------------------------------------------------------------ # 

G <- exp( G*Y_scale + Y_center ) # tranform to cm scale
R <- rstan::extract( rfit, 'IBM_Y_hat')$IBM_Y_hat

a <- get_recruit_area(spp)

K <- S*G  # survival by size 
R <- R*median(a)

K <- data.frame( quad = sd$IBM_quad_name, year = sd$IBM_year_name + 1, t(K))
R <- data.frame( quad = rd$IBM_quad_name, year = rd$IBM_year_name + 1, t(R))

K <- 
  K %>% 
  gather( sim, area, starts_with('X')) %>%
  group_by( quad, year, sim) %>% 
  summarise(area = sum( area ))

R <- 
  R %>% 
  gather( sim, area, starts_with('X')) 

A_pred <- 
  R %>% 
  left_join(K, by = c('quad', 'year', 'sim')) %>% 
  ungroup() %>% 
  gather( type, area, area.x, area.y) %>% 
  group_by( quad, year, sim) %>% 
  summarise( area = sum(area, na.rm = T)) %>%
  mutate( cover = 100*area/(100*100))

#---------------------------------------------------------- # 
# generate observed area per quad per year ---------------- # 
X <- 
  gdat %>% 
  mutate( quad =  QuadName) %>% 
  group_by( quad , year) %>% 
  summarise( area = sum(exp(logarea.t0), na.rm = T))  %>%
  mutate( cover0 = 100*area/(100*100))

K_obs <- gd$IBM_Y
K_obs <- exp( K_obs*Y_scale + Y_center)
R_obs <- median(a)*rd$IBM_Y

K_obs <- data.frame( quad = gd$IBM_quad_name, year = gd$IBM_year_name + 1, area = K_obs)
R_obs <- data.frame( quad = rd$IBM_quad_name, year = rd$IBM_year_name + 1, area = R_obs)

K_obs <- 
  K_obs %>% 
  group_by( quad, year ) %>% 
  summarise(area = sum( area, na.rm = T)) 

A_obs <- 
  R_obs %>% 
  left_join(K_obs, by = c('quad', 'year')) %>% 
  ungroup() %>% 
  gather( type, area, area.x, area.y) %>% 
  group_by( quad, year) %>% 
  summarise( area = sum(area, na.rm = T)) %>%
  mutate( cover = 100*area/(100*100))


A_obs <- 
  X %>% 
  full_join(A_obs , by = c('quad', 'year')) %>% 
  arrange( quad, year ) %>%
  ungroup( )%>% 
  mutate( cover_fill = ifelse( is.na(cover), cover0, cover)) %>%
  mutate( start_year = is.na(cover) & !is.na(cover0 )) %>% 
  select( quad, year, cover_fill, start_year)

A_obs <- 
  expand.grid( year = seq( min( A_obs$year) - 1, max(A_obs$year) ) + 1, quad = unique( A_obs$quad )) %>% 
  left_join(A_obs, by = c('quad', 'year')) 

A_pred <- 
  expand.grid( year = seq( min( A_pred$year) - 1, max(A_pred$year) ) + 1, quad = unique( A_pred$quad )) %>% 
  left_join(A_pred, by = c('quad', 'year')) %>% 
  left_join(quads, by = c('quad' = 'QuadName')) %>% 
  mutate( era = cut(year, include.lowest = T, breaks = c(min(year), 1960, 2004, 2018), labels = c('early', 'mid', 'late'))) %>% 
  select( quad, Treatment, year, era, cover, sim )

A_pred_summary <- 
  A_pred %>% 
  group_by( year, quad, Treatment, era) %>% 
  summarise( avg = mean(cover, na.rm = T), 
             low5 = quantile(cover, 0.05, na.rm = T), 
             low25 = quantile(cover, 0.25, na.rm = T),
             med50 = quantile(cover, 0.5, na.rm = T),
             upper75 = quantile(cover, 0.75, na.rm = T), 
             upper95 = quantile(cover, 0.95, na.rm = T)) %>% 
  mutate( year_label = as.numeric( str_sub(year, 3, 5)) )

pred_df <- 
  A_pred_summary %>% 
  left_join(
    A_obs %>% 
      rename('cover_obs' = cover_fill), 
    by = c('quad', 'year')
  )

saveRDS(pred_df, 'output/IBM_cover_predictions.RDS')

