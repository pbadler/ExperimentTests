library(lme4)
library(dplyr)
library(tidyr)

# call from removal_analysis_wrapper.r
root <- "~"
sppList=c("Artemisia tripartita","Hesperostipa comata","Poa secunda","Pseudoroegneria spicata")
dataDir <- file.path(root,"/driversdata/data/idaho_modern/")

load('analysis/figure_scripts/my_plotting_theme.Rdata')

# import data and calculate treatment trends ######################################

covD<-read.csv(paste(dataDir,"allrecords_cover.csv",sep=""))
trts<-read.csv(paste(dataDir,"quad_info.csv",sep=""))

# use this to make sure we don't miss zeros
allquadyrs<-unique(covD[,c("quad","year")],MARGIN=2)
tmp<-expand.grid(species=sppList,year=sort(unique(covD$year)))
allquadyrs<-merge(allquadyrs,tmp,all=T)
allquadyrs<-merge(allquadyrs,trts[,c("Group","quad","Treatment")])

# only do precip treatments
keep <- which(is.element(allquadyrs$Treatment,c("Control","Drought","Irrigation")))
allquadyrs <- allquadyrs[keep,]

# aggregate to quadrat level
keep<-which(is.element(covD$species,sppList))
sppD<-covD[keep,]
sppD<-merge(sppD,allquadyrs,all.y=T)
sppD$area[is.na(sppD$area)]<-0
sppD$area<-sppD$area*100
sppD.q<-aggregate(sppD$area,by=list(species=sppD$species,Group=sppD$Group,Treatment=sppD$Treatment,
                  quad=sppD$quad,year=sppD$year),FUN=sum)

test <- sppD.q %>% filter( species == 'Artemisia tripartita') %>%  group_by( species, year, Treatment ) %>% summarise( mean(x))

sample_size <- aggregate( sppD$area, by = list(species=sppD$species,Treatment=sppD$Treatment,
                                quad=sppD$quad,year=sppD$year), FUN = function(x) sum(x > 0))

sample_size_quad <- aggregate( sample_size$x, by = list(species=sample_size$species,Treatment=sample_size$Treatment,year=sample_size$year), FUN = function(x) sum(x > 0))
sample_size_count <- aggregate( sample_size$x, by = list(species=sample_size$species,Treatment=sample_size$Treatment,year=sample_size$year), FUN = sum)

sample_size_count$n <- sample_size_count$x
sample_size_quad$nquads <- sample_size_quad$x

sample_size_df <- merge( sample_size_quad, sample_size_count, by = c('species', 'Treatment', 'year'))

names(sppD.q)[NCOL(sppD.q)]<-"cover"
write.csv(sppD.q,"data/temp_data/QuadYearCover.csv",row.names = FALSE)  # save these 

sample_size_df <- sample_size_df[ , -grep('x', names(sample_size_df))]

sample_size_df <- sample_size_df[order(sample_size_df$species, sample_size_df$Treatment, sample_size_df$year), ]

write.csv(sample_size_df, 'output/sample_size_in_modern_data.csv')

# focus on all control plots or just those in the big exclosure?
#sppD.q <- subset(sppD.q,Group=="E1")
 
#calculate treatment means by year
spp.mean <- aggregate(sppD.q$cover,by=list(species=sppD.q$species,Treatment=sppD.q$Treatment,
                  year=sppD.q$year),FUN=mean)
names(spp.mean)[NCOL(spp.mean)] <- "cover"
spp.mean <- reshape(spp.mean,direction="wide",timevar="Treatment",idvar=c("species","year"))
spp.mean <- subset(spp.mean,year>2006)
spp.mean <- spp.mean[order(spp.mean$species,spp.mean$year),]

# calculate year-to-year log changes
tmp <- sppD.q
tmp$year <- tmp$year + 1
names(tmp)[which(names(tmp)=="cover")] <- "lag.cover"
logChange <- merge(sppD.q,tmp)
logChange$pcgr <- log(logChange$cover/logChange$lag.cover)
logChange$pcgr[logChange$pcgr==Inf] <- NA
logChange$pcgr[logChange$pcgr==-Inf] <- NA
mean.change <- aggregate(logChange$pcgr,by=list(species=logChange$species,Treatment=logChange$Treatment,
                  year=logChange$year),FUN=mean,na.rm=T)
names(mean.change)[NCOL(mean.change)] <- "pcgr"
mean.change <- subset(mean.change,year>2006)
mean.change <- reshape(mean.change,direction="wide",timevar="Treatment",idvar=c("species","year"))
mean.change <- mean.change[order(mean.change$species,mean.change$year),]

# calculate deviations from pretreatment year
sppD.q.2011<-subset(sppD.q,year==2011) #get pre-treatment year
names(sppD.q.2011)[NCOL(sppD.q.2011)]<-"cover.2011"
i<-which(names(sppD.q.2011)=="year")
sppD.q.2011<-sppD.q.2011[,-i]
sppD.q<-subset(sppD.q,year>2006)
sppD.q<-merge(sppD.q,sppD.q.2011)

sppD.q$coverDiff <- sppD.q$cover-sppD.q$cover.2011 # calculate difference 

spp.mean.diff<-aggregate(sppD.q$coverDiff,by=list(species=sppD.q$species,Treatment=sppD.q$Treatment,
                  year=sppD.q$year),FUN=mean)
names(spp.mean.diff)[NCOL(spp.mean.diff)]<-"coverDiff"
spp.mean.diff <- reshape(spp.mean.diff,direction="wide",timevar="Treatment",idvar=c("species","year"))
spp.mean.diff <- spp.mean.diff[order(spp.mean.diff$species,spp.mean.diff$year),]

# statistical tests ####################################################

dARTR <- subset(logChange,species=="Artemisia tripartita" & !is.na(pcgr) & Treatment!="No_shrub" & year > 2010 )
dARTR$year <- as.factor(dARTR$year)
dARTR$Treatment[ dARTR$year == 2011 & dARTR$Treatment != 'Control' ] <- 'Control'
mARTR <- lmer(pcgr ~ Treatment + (1|year),data=dARTR)

dHECO <- subset(logChange,species=="Hesperostipa comata" & !is.na(pcgr) & Treatment!="No_grass" & year > 2010 )
dHECO$year <- as.factor(dHECO$year)
dHECO$Treatment[ dHECO$year == 2011 & dHECO$Treatment != 'Control' ] <- 'Control'
mHECO <- lmer(pcgr ~ Treatment + (1|year),data=dHECO)

dPOSE <- subset(logChange,species=="Poa secunda" & !is.na(pcgr) & Treatment!="No_grass" & year > 2010 )
dPOSE$year <- as.factor(dPOSE$year)
dPOSE$Treatment[ dPOSE$year == 2011 & dPOSE$Treatment != 'Control' ] <- 'Control'
mPOSE <- lmer(pcgr ~ Treatment +  (1|year), data=dPOSE)

dPSSP <- subset(logChange,species=="Pseudoroegneria spicata" & !is.na(pcgr) & Treatment!="No_grass" & year > 2010 )
dPSSP$year <- as.factor(dPSSP$year)
dPSSP$Treatment[ dPSSP$year == 2011 & dPSSP$Treatment != 'Control' ] <- 'Control'
mPSSP <- lmer(pcgr ~ Treatment + (1|year),data=dPSSP)

species <- c('ARTR', 'HECO', 'POSE', 'PSSP')
statsOutput <- paste0( "output/", species, "_lmer_pgr_stats_table.text")
output <- list(mARTR,mHECO,mPOSE,mPSSP)

for(i in 1:length(species)){ 
  texreg::texreg(output[[i]], ci.force=TRUE,caption="Cover change models",
        caption.above=TRUE,file=statsOutput[i])
} 

summary(mARTR)
summary(mHECO)
summary(mPOSE)
summary(mPSSP)


# figures ########################################################################
library(ggplot2)
trtLabels<-substr(x=names(spp.mean)[3:5],start=7,stop=nchar(names(spp.mean)[3:5]))

trtLabels
spp.mean

spp.mean_cover_long <- 
  sppD.q %>% 
  filter( Group == 'E1') %>% 
  group_by( species, Treatment, year ) %>% 
  summarise( mcover = mean(cover))

spp.mean_cover_long$species <- factor( spp.mean_cover_long$species, labels = c('ARTR', 'HECO', 'POSE', 'PSSP'))

p1 <- 
  ggplot( spp.mean_cover_long, aes( x = year, y = mcover, color = Treatment ) ) + 
  geom_point() + 
  geom_line() + 
  ylab( 'Mean cover (%)' ) + 
  scale_color_manual(values = my_colors[2:4]) + 
  my_theme + theme(axis.text.x = element_text()) + 
  guides( color = 'none') + 
  scale_x_continuous(name = 'year', breaks = c(2007:2016)) 

p1 %+% subset(spp.mean_cover_long, species == 'ARTR') + ylim( 0, 25)

#1. Average cover treatment and year
png("figures/treatment_trends_cover.png",height=2.75,width=8,units="in",res=400)
  
  par(mfrow=c(1,4),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0),tcl=-0.2)

  # mean cover
  for(doSpp in sppList){
    tmp.mean<-subset(spp.mean,species==doSpp & year > 2010)
    if(doSpp==sppList[1]){
      my.y <- c(0,max(tmp.mean[,3:5], na.rm = T))
    }else{
      my.y<- c(0,3.2)
    }
    matplot(tmp.mean$year,tmp.mean[,3:5],ylim=my.y,type="o",xlab="",ylab="",pch=16,lty="solid",
            col=my_colors[2:4],main=doSpp,font.main=4,lwd=2)
    if(doSpp==sppList[1]) {
      legend("bottomright",c("Control","Drought","Irrigation"),pch=16,lty="solid",col=my_colors[2:4],bty="n")
      mtext("Mean cover (%)",side=2,line=2,outer=F)
    }
  }
  mtext("Year",side=1,line=1,outer=T)
  
dev.off()

#2. Log cover change by treatment and year 
png("figures/treatment_trends_logChange.png",height=3,width=8.5,units="in",res=400)
  
  par(mfrow=c(1,4),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2,0,0),tcl=-0.2,lwd=1)

  # log change
  myLims <- c(-2.2,2)
  for(doSpp in sppList){
    tmp.mean<-subset(mean.change,species==doSpp )
    #tmp.mean[1,3:5] <- NA  # get rid of control plot value in 2011
    #remove irrelevant removal treatments
    if(doSpp=="Artemisia tripartita"){
      tmp.mean$pcgr.No_shrub <- NA
    }else{
      tmp.mean$pcgr.No_grass <- NA
    }
    matplot(tmp.mean$year,tmp.mean[,3:5],type="o",xlab="",ylab="",pch=16,lty="solid",
            col=my_colors[2:4],ylim=myLims,main=doSpp, font.main=4)
    abline(h=0,col="gray")
    
    if(doSpp==sppList[1]) {
      legend("topright",c("Control","Drought","Irrigation"),pch=16,lty="solid",col=my_colors[2:4],bty="n")
      mtext("Mean log change",side=2,line=2,outer=F)
    }
    
  }
  mtext("Year",side=1,line=0.25,outer=T)

dev.off()

# 
#2. Average cover deviation (w.r.t. pretreatment year)
png("figures/cover_deviation.png", height=3,width=8.5,units="in",res=400)
par(mfrow=c(1,4),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2.5,0,0),tcl=-0.2)
for(doSpp in sppList){
  tmp.mean<-subset(spp.mean.diff,species==doSpp)
  matplot(tmp.mean$year,tmp.mean[,3:5],type="o",xlab="",ylab="",pch=16,lty="solid",
          col=my_colors[2:4],main=doSpp)

  if(doSpp==sppList[1]) legend("topleft",trtLabels,pch=16,lty="solid",col=my_colors[2:4],bty="n")
}
mtext("Year",side=1,line=1,outer=T)
mtext("Mean cover  deviation (%)",side=2,line=1,outer=T)
dev.off()
# 

# 3. Change in cover bar chart 

spp <- unique( sppD$species)

subset( sppD.q, species == 'Hesperostipa comata' & year == 2016) %>% 
  arrange( year ) %>% 
  filter( cover > 0 )


temp <- subset( sppD.q , year == 2016 )
temp$Treatment <- factor(temp$Treatment )
temp$lc <- log( temp$cover/temp$cover.2011 ) 
temp <- subset( temp , is.finite(lc))
temp$species <- factor( temp$species, labels =  c('ARTR', 'HECO', 'POSE', 'PSSP'))

lc.ARTR <- lm(dat =  subset( temp, species == 'ARTR'), lc ~ Treatment )
lc.HECO <- lm(dat =  subset( temp, species == 'HECO'), lc ~ Treatment )
lc.POSE <- lm(dat =  subset( temp, species == 'POSE'), lc ~ Treatment )
lc.PSSP <- lm(dat =  subset( temp, species == 'PSSP'), lc ~ Treatment )


summary(lc.ARTR)
summary(lc.HECO)
summary(lc.POSE)
summary(lc.PSSP)

library(xtable)
xt1 <- xtable(lc.ARTR, 
              caption = 'Treatment effects on log cover change for \textit{A. tripartita} from 2011 to 2016. Intercept gives control effects.', 
              label = 'table:changeARTR')
xt2 <- xtable(lc.HECO, caption = 'Treatment effects on log cover change for \textit{H. comata} from 2011 to 2016. Intercept gives control effects.', 
              label = 'table:changeHECO')
xt3 <- xtable(lc.POSE, caption = 'Treatment effects on log cover change for \textit{P. secunda} from 2011 to 2016. Intercept gives control effects.', 
              label = 'table:changePOSE')
xt4 <- xtable(lc.PSSP, caption = 'Treatment effects on log cover change for \textit{P. spicata} from 2011 to 2016. Intercept gives control effects.', 
              label = 'table:changePSSP')

print(xt1, 'output/results_tables/ARTR_cover_change.txt', type = 'latex')
print(xt2, 'output/results_tables/HECO_cover_change.txt', type = 'latex')
print(xt3, 'output/results_tables/POSE_cover_change.txt', type = 'latex')
print(xt4, 'output/results_tables/PSSP_cover_change.txt', type = 'latex')



temp$percent_change <- 100*(temp$coverDiff/temp$cover.2011 )

aggregate(data = temp, percent_change ~ Treatment + species, 'mean')

avgcover <- aggregate(data = temp, cbind(cover.2011, cover) ~ Treatment + species, 'mean')

avgcover$propchange <- (avgcover$cover - avgcover$cover.2011)/avgcover$cover.2011
avgcover

subset( temp , species == 'ARTR')

png("figures/start_to_finish_cover_change.png", height=4,width=10,units="in",res=400)
print( 
ggplot ( temp, aes( x  = Treatment, y = lc, color = Treatment )) + 
  geom_boxplot() + 
  geom_point( position = position_dodge(width = 1)) + 
  geom_hline(aes(yintercept = 0 ), linetype = 2) + 
  ylab( 'Log change in cover 2011 to 2016')  + 
  xlab ( '' ) + 
  facet_grid( . ~ species) + 
  scale_color_manual( values = my_colors[ 2:4]) + 
  my_theme
)
dev.off()

summary(m2[[1]])
summary(m2[[2]])
summary(m2[[3]])
summary(m2[[4]])

# #3. log change
# # hard wire ylims
# myLims <- c(-2.2,2)
# pdf("figures/log_change.pdf",height=3,width=8)
# par(mfrow=c(1,4),mgp=c(2,0.5,0),mar=c(2,2,2,1),oma=c(2,2.5,0,0),tcl=-0.2)
# for(doSpp in sppList){
#   tmp.mean<-subset(mean.change,species==doSpp)
#   #remove irrelevant removal treatments
#   if(doSpp=="Artemisia tripartita"){
#     tmp.mean$pcgr.No_shrub <- NA
#   }else{
#     tmp.mean$pcgr.No_grass <- NA
#   }
#   matplot(tmp.mean$year[2:5],tmp.mean[2:5,3:5],type="o",xlab="",ylab="",pch=16,lty="solid",
#           col=my_colors[2:4],main=doSpp,ylim=myLims)
#   abline(h=0,col="gray")
#   if(doSpp==sppList[1]) legend("top",trtLabels,pch=16,lty="solid",col=my_colors[2:4],bty="n")
# }
# mtext("Year",side=1,line=1,outer=T)
# mtext("Mean cover  deviation (%)",side=2,line=1,outer=T)
# dev.off()
