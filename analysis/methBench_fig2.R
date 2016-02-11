# FIGURE 2: COMPARISON OF ABSOLUTE ASSAYS


# ensure reproducibility by fixing the random seed to an arbitrary offset:
load(latestPrepName)
set.seed(1234)
message("=== FIGURE 2 ===")




### INTRODUCTION TO CORRIDOR CONCEPT (Figure 2A) ###

# plot an example beeswarm plots illustrating some corridors for region_04:
# (this an arbitrary example; all beeswarm plots are given as a supplementary figure!)
message("Example beeswarm panel for absolute assays (Figure 2.a)...")
plot3WayBeeswarm(lsaData[ lsaData$regionName=="region_04" & grepl("CRC_1|CRC_2|KG",lsaData$sampleName) & lsaData$datasetName%in%datasetsByType$absolute,],meth="swarm")




### PERFORMANCE FOR DIFFERENTIAL EXPERIMENTS (Figure 2B) ###

message("Boxplots for absolute and directional deviation from consensus corridor (Figure 2.b)...")
selDatasets <- datasetsByType$absolute
selRegions <- regionNames
selSamples <- setdiff(sampleNames,sampleNames[grep("Titration",sampleNames)])

# boxplots for absolute deviation to corridor:
message("\t* absolute deviation")
ggData <- melt(lsaData[lsaData$datasetName%in%selDatasets & lsaData$regionName%in%selRegions & lsaData$sampleName%in%selSamples,][,c("datasetName","deviationCorridor")],id.vars="datasetName")
ggData <- cbind(ggData,"col"=makeTransparent(plotColLookup$datasetName[ggData$datasetName],0.75))
ggData$datasetName <- factor(datasetTable[ggData$datasetName,"prettyLabel"],levels=rev(datasetTable[datasetNamesInOrder,"prettyLabel"]))
ggData$variable <- "Relative to consensus corridor"
d <- ggplot(ggData, aes(x=datasetName, y=abs(value), fill=col)) + geom_boxplot(notch=FALSE,outlier.size=1.5,type="count") +ylim(0,100) + scale_fill_identity(guide = FALSE) + xlab(NULL) +coord_flip()+ ylab("Absolute deviation") + defaultPlotTheme(flipX=TRUE)
svgPlotGG(d,"performance_differential", 9, 5.8, units="in")

# boxplots for "bias" (i.e. trend to over- or under-estimate value):
message("\t* directional deviation")
ggData2 <- melt(lsaData[lsaData$datasetName%in%selDatasets & lsaData$regionName%in%selRegions & lsaData$sampleName%in%selSamples,][,c("datasetName","deviationCorridor")],id.vars="datasetName")
ggData2 <- aggregate(value ~ datasetName + variable, ggData2, mean, simplify=T)
ggData2 <- cbind(ggData2,"col"=makeTransparent(plotColLookup$datasetName[ggData2$datasetName],0.75))
ggData2$datasetName <- factor(datasetTable[ggData2$datasetName,"prettyLabel"],levels=rev(datasetTable[datasetNamesInOrder,"prettyLabel"]))
ggData$variable <- "Relative to consensus corridor"
d <- ggplot(ggData2, aes(x=datasetName,y=value,fill=value)) + geom_bar(stat="identity",colour="#333333") +coord_flip() + ylim(-5,5) + xlab(NULL) + ylab("Bias")  + scale_fill_gradient2(low=brewer.pal(3,"BrBG")[3],high=brewer.pal(3,"BrBG")[1],guide=F,space="Lab") + geom_hline(aes(yintercept=0),colour="#333333") + defaultPlotTheme(flipX=TRUE)  
svgPlotGG(d,"performance_differential_bias_quant", 2.8, 5.8, units="in")



### INTRODUCTION TO TITRATION SERIES (Figure 2C) ###

message("Example plots for titration series (Figure 2.c)...")

# pick the best- and worst-performing assays:
# (selection criteria: only assays that attempted both titration series, amongst those: the highest/lowest *median* Pearson correlation coefficient)
exclDatasets <- c("Pyroseq_5","Pyroseq_1_replicate","AmpliconBS_3") # these assays have only attempted one titration series, so we exclude them here
tmp <- performanceTableTitrations[performanceTableTitrations$datasetName %in% setdiff(datasetsByType$absolute,exclDatasets),]; colnames(tmp)[5] <- "Pearson"; tmp <- tmp[!is.na(tmp$Pearson),];tmp <- aggregate(Pearson~datasetName,tmp,median); tmp <- tmp[order(tmp$Pearson,decreasing=T),]
bestTitrationAssay <- tmp[1,"datasetName"]
worstTitrationAssay <- tmp[nrow(tmp),"datasetName"]
rm(tmp)

selDatasets <- c(bestTitrationAssay,worstTitrationAssay) 
selTitrations <- c("Titration1","Titration2")
svgPlot("titration_analysis_panel", 6, 6, pointsize=14)
par(mfrow=c(2,2))
par(mgp=c(1.75,0.5,0))
par(mar=c(3,3,2,1)+0.1)
for (titrationType in selTitrations) {
	message("\t* ", titrationType)
	curTitrationData <- correctedTitrationData[correctedTitrationData$titrationType==titrationType,]
	rng <- c(-25,100) # fix plot range

	for (curDataset in selDatasets) {
		curData <- curTitrationData[curTitrationData$datasetName==curDataset,]
		
		isLeftCol <- curDataset==selDatasets[1]
		isTopRow <- titrationType==selTitrations[1]
		par(mar=c(ifelse(isTopRow,1,4),ifelse(isLeftCol,3,1),ifelse(isTopRow,4,1),ifelse(isLeftCol,1,3))+0.1)

		curData$x <- curData$titrationPercent 
		# perform log-transformation for better visualization (for titration series 2):
		doLog <- titrationType=="Titration2"
		if(doLog) { curData$x <- log10(curData$x*1000+1) }
		xmax <- max(curData$x,na.rm=T)

		# fit a model to the offset-adjusted data:
		fit <- lm(methAdj ~ titrationPercent - 1,data=curData) # no offset, it's already subtracted!
		curData$prediction <- fitted(fit)

		plot(curData$x,curData$methAdj, xlim=c(min(curData$x,na.rm=T),xmax), ylim=rng, xlab="Titration percentage", las=1, ylab=ifelse(isLeftCol,"Adj. DNA methylation (%)",""), main=NA, bty="l", col="black", cex=0.5, xaxt=ifelse(doLog,"n","s"))
		
		# if log-adjusted, add a corresponding x-axis:
		if(doLog) {
			ticks <- ((10^curData$x)-1)/1000
			labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
			axis(1, at=(1:length(ticks))-1, labels=rev(ticks))
		}

		# add a line per region:
		for (regionName in selRegions) {
			sel <- curData$regionName==regionName
			lines(curData$x[sel],curData[sel,"methAdj"],lwd=0.5)
		}

		# add a line for the fitted model:
		if (length(coef(fit))<=2) { 
			x <- 0:100
			xx <- x
			if(doLog) xx <- log10(x*1000+1)
			lines(xx,predict(fit,data.frame("titrationPercent"=x)),col="hotpink",lwd=3)
		}
		else cat("No fit available?\n")

		# add labels:
		if(!isLeftCol) {
			text(xmax*1.12,sum(rng)/2,titrationType,srt=270,xpd=NA,cex=1.5,font=2)
		}
		if(isTopRow) {
			text(50,rng[2]*1.4,curDataset,xpd=NA,cex=1.5,font=2)
		} 

		# calculate the median Pearson correlation between the measurements for a region and the expected target methylation levels:
		r <- median(performanceTableTitrations[performanceTableTitrations$datasetName==curDataset & performanceTableTitrations$titrationType==titrationType,"Pearson.s.r"])
		text(0,rng[2],paste("r =",format(r,digits=3,nsmall=3)),adj=c(0,1))
		
	}
}
dev.off()




### PERFORMANCE FOR TITRATION SERIES (Figure 2D) ###

message("Performance characteristics for titration series (Figure 2.d)...")

# prepare data by melting the data table and adjusting labels:
ggData <- performanceTableTitrations[
	regexpr("region_",performanceTableTitrations$regionName)>=0 & performanceTableTitrations$datasetName%in%datasetsByType$absolute,
	c("datasetName","titrationType","Pearson.s.r","Adjusted.R.2","Residual.standard.error")
]
ggData <- melt(ggData,id.vars=c("datasetName","titrationType"))
ggData <- cbind(ggData,"col"=makeTransparent(plotColLookup$datasetName[ggData$datasetName],0.75))
ggData$datasetName <- datasetTable[ggData$datasetName,"prettyLabel"]
levels(ggData$variable)[levels(ggData$variable)=="Pearson.s.r"] <- "Pearson's r"
levels(ggData$variable)[levels(ggData$variable)=="Adjusted.R.2"] <- "Adjusted R^2"
levels(ggData$variable)[levels(ggData$variable)=="Residual.standard.error"] <- "Residual standard error"

# plot a panel of boxplots:
d <- ggplot(ggData, aes(datasetName,value,fill=col)) + geom_boxplot(outlier.size=1) + facet_wrap(titrationType ~ variable, scales="free_y") + xlab("") + ylab("") + scale_fill_identity(guide = FALSE) + defaultPlotTheme(flipX=TRUE)
# ... add frame for all plots in panel:
d <- d +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
svgPlotGG(d,"performance_titrations", 7, 4.6, units="in")





save.image(paste("methBench","afterFig2",format(Sys.time(), "%Y%m%d_%H%M"),"session.Rdata",sep="_"))

