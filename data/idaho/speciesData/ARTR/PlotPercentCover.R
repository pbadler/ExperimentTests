setwd("c:/repos/driversData/data/idaho/speciesData/ARTR");
X = read.csv(file="quadratCover.csv"); 

years=unique(X$year); 
meanCover = numeric(length(years))  
for(j in 1:length(years)) {
	meanCover[j] = median(X$totCover[X$year==years[j]])
}

pctCover = meanCover/(100*100); 
plot(years,pctCover,type="b",xlab="Year",ylab="Median %cover of sampled quadrats"); 

xvals=years[-(1:6)]; yvals=mean(pctCover[-(1:6)]); yvals=rep(yvals,length(xvals));
points(xvals,yvals,lty=2,type="l"); 
title(main=paste("ARTR: average of annual median covers (w/o 1st 6 years)=",round(100*yvals[1],digits=1))); 