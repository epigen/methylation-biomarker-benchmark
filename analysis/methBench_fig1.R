# FIGURE 1: SCATTERPLOTS AND CLUSTERING


# ensure reproducibility by fixing the random seed to an arbitrary offset:
load(latestPrepName)
set.seed(1234)
message("=== FIGURE 1 ===")




### CORRELATION HEATMAPS (Figure 1C) ###

# perform hierarchical clustering and plot heatmap:
message("Correlation heatmaps (Figure 1.c)...")
for (curType in c("absolute","relative")) { 
	curDatasets <- datasetsByType[[curType]]

	message("\t* ",curType)
	curCorrelations <- corTable[curDatasets,curDatasets]
	colnames(curCorrelations) <- datasetTable[curDatasets,"prettyLabel"]
	rownames(curCorrelations) <- datasetTable[curDatasets,"prettyLabel"]
	clust <- hclust(dist(1-abs(curCorrelations)),"average")

	svgPlot(paste0("cor_hierarchical_clustering_",curType), 5, 5)
	par(mar=c(4,2,1,8))
	dend <-as.dendrogram(clust)
	plot(dend,horiz=TRUE)
	dev.off()
	
	cOrder <- rank(clust$order[datasetNames %in% curDatasets])

	svgPlot(paste("cor_heatmap_",curType), 12, 12)
	pheatmap(curCorrelations[cOrder,cOrder],col=colorRampPalette(rev(brewer.pal(5,"Oranges")),bias=0.2,space="Lab")(20),breaks=seq(0,1,0.05),border_color="white",cellwidth=22,cellheight=22,cluster_cols=FALSE,cluster_rows=FALSE,display_numbers=TRUE,fontsize=13,fontcolor="black")
	dev.off()
}




### EXAMPLE SCATTERPLOTS (Figure 1D) ###

# for selection of datasets: 
# sort(corTable[rownames(corTable) %in% datasetsByType$absolute,"Pyroseq_1"])

message("Example scatter plots (Figure 1.d)...")
for(s2 in c("Pyroseq_1")) {
	for(s1 in c("Pyroseq_1_replicate","Pyroseq_2","Infinium","EpiTyper_3","AmpliconBS_1","EnrichmentBS_1")) {
		message("\t* ",s1,"-vs-",s2)
		ggDataSelected <- assayVsAssayData[assayVsAssayData$assay1==s1&assayVsAssayData$assay2==s2,]
		d <- ggplot(data=ggDataSelected,aes(x=meth1,y=meth2)) + geom_point(shape=16,alpha=0.3) + geom_smooth(method=lm,fill='steelblue') + theme_bw() + xlab(DNA_METH_LABEL) + ylab(DNA_METH_LABEL) + xlim(0,100) + ylim(0,100) + annotate("text", x = 100, y = 0, label = paste("r =",format(cor(ggDataSelected$meth1,ggDataSelected$meth2),digits=3,nsmall=3)),hjust=1,vjust=0,size=4) + labs(title = datasetTable[s1,"prettyLabel"]) + theme(			axis.line = element_line(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank()
		)		 
		
		svgPlotGG(d, paste0("individual_scatterplots_ggplot_",s1,"-vs-",s2), 3, 3, units="in")
	}
}




### BEESWARM OVERVIEW PLOTS (Figure 1B) ###

message("Beeswarm overview plots (Figure 1.b)...")
svgPlot("beeswarm_1D_panel", 5, 3.4, pointsize=11)
 
temp <- lsaData[lsaData$datasetName%in%datasetsByType[["absolute"]],]
temp <- temp[sample(nrow(temp),nrow(temp)),] # randomize order so not to bias the looks of the plots

par(xpd=NA)
par(mfrow=c(3,1))
par(mgp=c(1.75,0.5,0))
par(mar=c(0.25,2.75,2,0.5)+0.1)

# plots the three components of the figure panel:
for(curCat in c("dataset","sample","region")) {
	message("\t* ",curCat)
	bees <- beeswarm(
			as.formula(paste("methValue","~",curCat,"Name",sep="")), data=temp, labels=NA,
			horizontal=FALSE, cex=0.05, spacing=2, cex.lab=1, pch=16, col="#555555",method="hex",corral="random",main=paste(length(unique(temp[,paste(curCat,"Name",sep="")]))," ",curCat,"s",sep=""),xlab=NA,ylab=NA,ylim=c(0,100),pwcol=makeTransparent(plotColLookup$datasetName[temp$datasetName],0.75), bty='l', xaxt='n',yaxt='n'
	)
	axis(2,seq(0,100,25),as.character(c(0,NA,50,NA,100)),las=1,cex=0.9)
}

dev.off() 

# p-values (for pairwise differences in overall distributions):

temp2 <- temp[order(temp$datasetName,temp$regionName,temp$sampleName),]
cmps <- list(
	"TumorNormal"=c("Tumor","Normal","two.sided"),
	"DrugControl"=c("Control","AzaC","two.sided"),
	"FrozenFFPE"=c("FFPE","Frozen","two.sided"),
	"Titration1"=c("1_100","6_0","two.sided"),
	"Titration2"=c("1_100","6_0","two.sided")
)
pvals <- c()
for(selSampleType in names(cmps)) {
	g1 <- grep(cmps[[selSampleType]][1],sampleNamesByType[[selSampleType]],value=TRUE)
	g2 <- grep(cmps[[selSampleType]][2],sampleNamesByType[[selSampleType]],value=TRUE)
	temp3 <- temp2[temp2$sampleName%in%c(g1,g2),]
	temp3$simpleSampleName <- gsub(paste0("(",cmps[[selSampleType]][1],"|",cmps[[selSampleType]][2],")"),"",temp3$sampleName)
	temp3g1 <- temp3[temp3$sampleName%in%g1,c("simpleSampleName","regionName","datasetName","methValue")]
	temp3g2 <- temp3[temp3$sampleName%in%g2,c("simpleSampleName","regionName","datasetName","methValue")]
	temp3mrg <- merge(temp3g1,temp3g2, by=c("simpleSampleName","regionName","datasetName"), all=FALSE)
	pval <- t.test(beta2mval(temp3mrg$methValue.x/100), beta2mval(temp3mrg$methValue.y/100), alternative =  cmps[[selSampleType]][3], paired = TRUE)
	message("paired, two-sided t-test for ", selSampleType, ": ", cmps[[selSampleType]][1], " vs ", cmps[[selSampleType]][2], " (",cmps[[selSampleType]][3],")?\tp-value = ", pval$p.value)
	pvals <- c(pvals,pval$p.value)
	names(pvals)[length(pvals)] <- selSampleType
}
print(round(p.adjust(pvals),8))

# same with corrected titration data:

temp2 <- correctedTitrationData[order(correctedTitrationData$datasetName,correctedTitrationData$regionName,correctedTitrationData$titrationType,correctedTitrationData$titrationPercent),]
temp2 <- temp2[temp2$titrationPercent%in%c(0,100),]
cmps <- list(
	"Titration1"=c(100,0,"two.sided"),
	"Titration2"=c(100,0,"two.sided")
)
pvals <- c()
for(titrationType in names(cmps)) {
	g1 <- temp2$titrationPercent == cmps[[titrationType]][1] & temp2$titrationType == titrationType
	g2 <- temp2$titrationPercent == cmps[[titrationType]][2] & temp2$titrationType == titrationType
	temp3g1 <- temp2[g1,c("regionName","datasetName","meth")]
	temp3g2 <- temp2[g2,c("regionName","datasetName","meth")]
	temp3mrg <- merge(temp3g1,temp3g2, by=c("regionName","datasetName"), all=FALSE)
	pval <- t.test(beta2mval(temp3mrg$meth.x/100), beta2mval(temp3mrg$meth.y/100), alternative = cmps[[titrationType]][3], paired = TRUE)
	message("paired, two-sided t-test for ", titrationType, ": ", cmps[[titrationType]][1], " vs ", cmps[[titrationType]][2], " (",cmps[[titrationType]][3],")?\tp-value = ", pval$p.value)
	pvals <- c(pvals,pval$p.value)
	names(pvals)[length(pvals)] <- titrationType
}
print(round(p.adjust(pvals),8))




save.image(paste("methBench","afterFig1",format(Sys.time(), "%Y%m%d_%H%M"),"session.Rdata",sep="_"))
