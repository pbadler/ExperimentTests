# format prediction tables 
alllppd <- read.csv('output/overall_lppd.csv')
allMSE <- read.csv('output/all_MSE_scores.csv')

allMSE$score <- 'MSE'
alllppd$score <- 'lppd'

alllppd$model <- factor(alllppd$model, labels = c('climate model', 'year effects model'))
allMSE$model <- factor(allMSE$model, labels = c('climate model', 'year effects model'))

all_scores <- rbind( allMSE %>% spread(model, MSE), alllppd %>% spread(model, lppd) )
all_scores$diff <- all_scores$`climate model` - all_scores$`year effects model`

overall_score_table <- 
  all_scores %>% 
  mutate( improved = ifelse( diff < 0 & score == 'MSE', '***', '')) %>% 
  mutate( improved = ifelse( diff > 0 & score == 'lppd', '***', improved)) %>% 
  arrange( species, vital_rate, score )


oaxt <- xtable(overall_score_table, caption = 'Comparison of model predictions from climate model and year effects model for each species and vital rate.  Two prediction scores are reported, MSE and lppd. Lower MSE indicates improved predictions whereas higher lppd indicates improved predictions.  Instances where the climate model outperformed the random year effects model are marked with "***" in the last column. ARTR = \\textit{A. tripartita}, HECO = \\textit{H. comata}, POSE = \\textit{P. secunda}, PSSP = \\textit{P. spicata}.' , 
label = 'table:overallPreds')

print(oaxt, 'manuscript/overall_prediction_score.tex', type = 'latex', include.rownames = F, caption.placement ="top")

# ---------- predictions by treatment level ------------------------------------------------------- # 

tlppd <- read.csv('output/treatment_lppd.csv')
tMSE <- read.csv('output/treatment_MSE_scores.csv')

tlppd <- tlppd %>% dplyr::select(-climate_model_performance) %>% gather( model, lppd, climate:year)

tMSE$score <- 'MSE'
tlppd$score <- 'lppd'


tlppd$model <- factor(tlppd$model, labels = c('climate model', 'year effects model'))
tMSE$model <- factor(tMSE$model, labels = c('climate model', 'year effects model'))


treatment_scores <- rbind( tMSE %>% spread(model, MSE), tlppd %>% spread(model, lppd) )
treatment_scores$diff <- treatment_scores$`climate model` - treatment_scores$`year effects model`

treatment_score_table <- 
  treatment_scores %>% 
  mutate( improved = ifelse( diff < 0 & score == 'MSE', '***', '')) %>% 
  mutate( improved = ifelse( diff > 0 & score == 'lppd', '***', improved)) %>% 
  arrange( species, vital_rate, Treatment, score )

treatment_score_table <- rename( treatment_score_table , `no climate model` = `year effects model`)

txt <- xtable(treatment_score_table, caption = 'Comparison of model predictions from climate model and no climate model for each species, vital rate and treatment.  Two prediction scores are reported, MSE and lppd. Lower MSE indicates improved predictions whereas higher lppd indicates improved predictions.  Instances where the climate model outperformed the no climate model are marked with "***" in the last column. ARTR = \\textit{A. tripartita}, HECO = \\textit{H. comata}, POSE = \\textit{P. secunda}, PSSP = \\textit{P. spicata}.' , 
               label = 'table:treatmentPreds')


add.to.row <- list(pos = list(0), command = NULL)
command <- paste0("\\hline\n\\endhead\n",
                  "\\hline\n",
                  "\\multicolumn{", dim(txt)[2] + 1, "}{l}",
                  "{\\footnotesize Continued on next page}\n",
                  "\\endfoot\n",
                  "\\endlastfoot\n")
add.to.row$command <- command

print(txt, 'manuscript/treatment_prediction_score.tex', type = 'latex', show.rownames = F, hline.after=c(-1), add.to.row = add.to.row,
      tabular.environment = "longtable", floating = F, caption.placement ="top") 


species_improvement <- 
  treatment_score_table %>% 
  group_by( species, Treatment, score ) %>% 
  summarise( sum( improved == '***')/ n() )

species_improvement

treatment_improvement <- 
  treatment_score_table %>% 
  group_by( Treatment, score ) %>% 
  summarise( improved_n =  sum( improved == '***'), total = n() )

treatment_improvement


treatment_vr_improvement <- 
  treatment_score_table %>% 
  group_by( Treatment, vital_rate, score ) %>% 
  summarise( improved_n =  sum( improved == '***'), total = n() )

treatment_vr_improvement
