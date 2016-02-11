# FIGURE 3: COMPARISON OF RELATIVE ASSAYS


# ensure reproducibility by fixing the random seed to an arbitrary offset:
load(latestPrepName)
set.seed(1234)
message("=== FIGURE 3 ===")



# define matched comparisons:
comparisons <- data.frame(
	"FG"=c(sampleNames[seq(2,12,by=2)],sampleNames[seq(13,16,by=2)],sampleNames[seq(29,32,by=2)]),
	"BG"=c(sampleNames[seq(1,12,by=2)],sampleNames[seq(14,16,by=2)],sampleNames[seq(30,32,by=2)])
,stringsAsFactors=FALSE)

# define bins and thresholds for discretizing consensus differences:
BINS <- list(
	c(0,5),
	c(5,25),
	c(25,50),
	c(50,101)
)
MARGINALITY_THRESHOLD_RELATIVE <- 0  # only absolute differences of 0 are to be considered as "no change observed" for relative assays
MARGINALITY_THRESHOLD_ABSOLUTE <- BINS[[1]][2] # for absolute assays, all changes smaller than a threshold of 5 are "no change observed"
names(BINS) <- sapply(BINS,function(x) { paste(x[1],"-",x[2],sep="")})
names(BINS)[1] <- paste("<",BINS[[1]][2])
names(BINS)[length(BINS)] <- paste(">",BINS[[length(BINS)]][1])

# select data to work with:
selRegions <- regionNames
selSamples <- setdiff(sampleNames,sampleNames[grep("Titration",sampleNames)])
selDatasets <- datasetNames
temp1 <- lsaData[lsaData$datasetName%in%selDatasets & lsaData$regionName%in%selRegions & lsaData$sampleName%in%selSamples,]

# calculate differences between consensus values and compare these with the changes recorded by the relative assays:
message("Perform comparisons between relative assays and consensus difference...")
performanceTable <- numeric()
for(cmpI in 1:nrow(comparisons)) {
	cmp <- as.character(comparisons[cmpI,])
	cat("\t*",paste(cmp,collapse=" vs "),"\n")

	# get consensus values for fore- and background:
	consensusFG <- consensus$median[selRegions,cmp[1]]
	consensusBG <- consensus$median[selRegions,cmp[2]]

	# calculate absolute differences:
	consensusMedianDiff <- consensus$median[selRegions,cmp[1]]-consensus$median[selRegions,cmp[2]]
	# ... set the difference to 0 if the consensus corridors overlap:
	consensusDiff <- ifelse(consensus$upper[selRegions,cmp[1]]>=consensus$lower[selRegions,cmp[2]] & consensus$lower[selRegions,cmp[1]]<=consensus$upper[selRegions,cmp[2]], 0, consensusMedianDiff)

	# find the bins the consensus differences fits into (should be only one, since they are non-overlapping):
	fittingBins <- sapply(BINS,function(x) { abs(consensusDiff)<x[2] & abs(consensusDiff)>=x[1]}, simplify=T)
	
	# then proceed one region at a time...
	for(curRegion in selRegions) {
		# get the data for the current region:
		temp2 <- temp1[temp1$regionName==curRegion & temp1$sampleName %in% cmp, c("datasetName","sampleName","methValue")]
		if(nrow(temp2)==0) next
		curConsMedDiff <- consensusMedianDiff[curRegion]
		curConsDiff <- consensusDiff[curRegion]
		curConsDir <- sign(curConsMedDiff)
		if(is.na(curConsDiff)) next
		
		# check for each dataset how it performs:
		for(curDataset in selDatasets) {		
			# get the methylation values for the current dataset:
			temp3 <- temp2[temp2$datasetName==curDataset, c("sampleName","methValue")]
			relFG <- temp3[temp3$sampleName==cmp[1],"methValue"]
			relBG <- temp3[temp3$sampleName==cmp[2],"methValue"]
			if(length(relFG)+length(relBG)<2) next
				
			# calculate difference:
			relDiff <- relFG - relBG

			# define "marginality" depending on whether it's an absolute or a relative assay that we're looking at:
			MARGINALITY_THRESHOLD <- ifelse(curDataset%in%datasetsByType$absolute,MARGINALITY_THRESHOLD_ABSOLUTE,MARGINALITY_THRESHOLD_RELATIVE)

			# make a call: either "no change detected" or the "same"/"opposite" as the consensus:
			ok <- ifelse(abs(relDiff)<=MARGINALITY_THRESHOLD,"Undetected",ifelse(sign(relDiff)==curConsDir,"Same","Opposite"))
			
			# remember performance metrics:
			performanceTable <- rbind(performanceTable,data.frame(
				"dataset" = curDataset,
				"region" = curRegion,	 			
				"titrationType" = paste(cmp,collapse=" x "),
				"expDifference" = which(fittingBins[curRegion,]),
				"call" = ok,
				"relativeChange" = ifelse(abs(relDiff)<=MARGINALITY_THRESHOLD,0,sign(relDiff)),				
				"magnitudeOfChange" = abs(relDiff)
			))
		}	
	}
}	
relativePerformanceTable <- performanceTable




### ASSESSMENT OF ABILITY OF RELATIVE ASSAYS TO DETECT CORRECT DIRECTION OF CHANGE (Figure 3C (first part: differential comparisons)) ###

# relative to consensus corridor differences:
message("Relative assay vs. consensus change plots (Figure 3.c)...")

# prepare data:
ggData <- performanceTable[performanceTable$dataset %in% datasetsByType$relative,]
ggData$dataset <- datasetTable[as.character(ggData$dataset),"prettyLabel"]
ggData$expDifference <- factor(names(BINS)[ggData$expDifference],names(BINS))
ggData$call <- as.character(ggData$call)

# how many of each type
smallChangeCalls <- aggregate(relativeChange~call+dataset,ggData[ggData$expDifference=="5-25",],function(x) { length(x)  })
aggregate(relativeChange~dataset,smallChangeCalls,function(x) { x[1]/sum(x) })

d <- ggplot(ggData, aes(x=expDifference,fill=call)) + geom_bar(width=0.5,colour="#333333",size=0.25)  + scale_fill_manual(values=c("Opposite"="firebrick","Same"="steelblue","Undetected"="gray")) + facet_grid(. ~ dataset, scales="fixed", margins=F) + coord_flip()   + xlab("Consensus difference")  + ylab("Number of observations") + geom_hline(yintercept=0) + guides(fill = guide_legend(title=NULL,keywidth=0.5,keyheight=0.5,label.theme = element_text(size = 8,angle=0))) + defaultPlotTheme()+ theme(axis.line.x = element_line(size=0.22), axis.line.y=element_blank())
svgPlotGG(d, "performance_rel_v_abs", 18, 4, units="in") 




### ASSESSMENT OF HOW "QUANTITATIVE" RELATIVE ASSAYS ARE (Figure 3D) ###

message("Boxplots measuring quantitative-ness (Figure 3.d)...")
ggData <- performanceTable[performanceTable$dataset %in% c("Pyroseq_1",datasetsByType$relative),]
ggData$dataset <- datasetTable[as.character(ggData$dataset),"prettyLabel"]
ggData$expDifference <- factor(names(BINS)[ggData$expDifference],names(BINS))
ggData$call <- as.character(ggData$call)
ggData2 <- ggData[ggData$call=="Same",] # only use the measurements that went into the correct direction of change!

d <- ggplot(ggData2, aes(dataset,magnitudeOfChange,fill=expDifference)) + xlab(NULL) + ylab("Measured difference") + guides(fill = guide_legend(title="Consensus\ndifference",keywidth=0.9,keyheight=0.9,label.theme = element_text(size = 8,angle=0))) + scale_fill_brewer() + scale_x_discrete(limits=rev(unique(ggData2$dataset))) + geom_boxplot(outlier.size=1) + coord_flip() + defaultPlotTheme()
svgPlotGG(d, "performance_rel_v_abs_quantitativity2", 11, 8)




### ASSESSMENT OF ABILITY OF RELATIVE ASSAYS TO DETECT CORRECT DIRECTION OF CHANGE (Figure 3C (second part: titration series)) ###

# relative to know direction of change in titrations series:
message("Perform comparisons between relative assays and titration series...")

selTitrations <- c("Titration1","Titration2")
selRegions <- regionNames
selDatasets <- datasetNames
performanceTable <- numeric()

for(titrationType in selTitrations) {
	cmp <- as.character(comparisons[cmpI,])
	cat("\t*",titrationType,"\n")
	
	# get the data for the current titration:
	titrationData <- correctedTitrationData[correctedTitrationData$titrationType==titrationType & correctedTitrationData$datasetName %in% selDatasets,]
	titrationData <- split(titrationData,titrationData$titrationPercent)
	
	# get the titration levels, sorted:
	titrationLevels <- as.character(sort(as.numeric(names(titrationData))))
	
	BASELINE_INDEX <- 1

	for(titrationLevelI in 1:length(titrationLevels)) {
		if(titrationLevelI == BASELINE_INDEX) next;
		
		# calculate the expected difference in DNA methylation as the differences in titration percentage:
		FG <- titrationLevels[titrationLevelI]
		BG <- titrationLevels[BASELINE_INDEX]
		expectedDiff <- as.numeric(FG) - as.numeric(BG)
		expectedDir <- sign(expectedDiff)

		# get the actual measurements:
		curTitrationLevel <- titrationData[[FG]]
		prevTitrationLevel <- titrationData[[BG]]

		# define "marginality" depending on whether it's an absolute or a relative assay that we're looking at:
		MARGINALITY_THRESHOLD <- ifelse(curTitrationLevel$datasetName%in%datasetsByType$absolute,MARGINALITY_THRESHOLD_ABSOLUTE,MARGINALITY_THRESHOLD_RELATIVE)

		# calculate observed differences:
		diffs <- curTitrationLevel$meth - prevTitrationLevel$meth

		# make a call: either "no change detected" or the "same"/"opposite" as the consensus:
		ok <- ifelse(abs(diffs)<=MARGINALITY_THRESHOLD,"Undetected",ifelse(sign(diffs)==expectedDir,"Same","Opposite")) 

		
		fittingBins <- sapply(BINS,function(x) { abs(expectedDiff)<x[2] & abs(expectedDiff)>=x[1]}, simplify=T)

		# remember performance metrics:
		performanceTable <- rbind(performanceTable,data.frame(
			"dataset" = curTitrationLevel$datasetName,
			"region" = curTitrationLevel$regionName,	 
			"titrationType" = titrationType,
			"expDifference" = as.numeric(which(fittingBins)),
			"call" = ok,
			"relativeChange" = ifelse(abs(diffs)<=MARGINALITY_THRESHOLD,0,sign(diffs)),
			"magnitudeOfChange" = abs(diffs)
		))

		# also remember the performance metrics "per titration series":
		performanceTable <- rbind(performanceTable,data.frame(
			"dataset" = curTitrationLevel$datasetName,
			"region" = curTitrationLevel$regionName,	 
			"titrationType" = paste(ifelse(titrationType=="Titration1","T1","T2"),": ",as.numeric(FG),"<>",as.numeric(BG),collapse=""),
			"expDifference" = as.numeric(which(fittingBins)),
			"call" = ok,
			"relativeChange" = ifelse(abs(diffs)<=MARGINALITY_THRESHOLD,0,sign(diffs)),
			"magnitudeOfChange" = abs(diffs)
		))
	}
}	
relativePerformanceTable <- rbind(relativePerformanceTable,performanceTable)


# plot the results of the titration series analysis:

message("Relative assay vs. consensus change plots (Figure 3.c, second part: titrations)...")

ggData <- performanceTable[performanceTable$dataset %in% datasetsByType$relative,]
ggData$dataset <- datasetTable[as.character(ggData$dataset),"prettyLabel"]
ggData$call <- as.character(ggData$call)
ggData$titrationType <- factor(ggData$titrationType,levels=rev(unique(ggData$titrationType)))

d <- ggplot(ggData, aes(x=titrationType,fill=call)) + geom_bar(width=0.5,colour="#333333",size=0.25)  + scale_fill_manual(values=c("Opposite"="firebrick","Same"="steelblue","Undetected"="gray")) + facet_grid(. ~ dataset, scales="fixed", margins=F) + coord_flip()   + xlab("Consensus difference\n")  + ylab("Number of observations") + geom_hline(yintercept=0) + guides(fill = guide_legend(title=NULL,keywidth=0.5,keyheight=0.5,label.theme = element_text(size = 8,angle=0))) + defaultPlotTheme(fontSize=8) + theme(axis.line.x = element_line(size=0.22), axis.line.y=element_blank())
svgPlotGG(d, "performance_rel_v_titrations", 18, 7)



### RELATIVE ASSAY COMPARISON (Figures 3A+B) ###


message("Comparisons of relative assays (Figures 3.a and 3.b)...")

message("\t* calculations")
# select data:
selRegions <- unique(relativePerformanceTable$region)
tTypes <- unique(relativePerformanceTable$titrationType)
selDatasets <- datasetNames
selSamples <- tTypes[1:10] # focus on the differences in the pair-wise comparisons, which have been attempted by all assays and give larger, more robust differences that make for a fairer comparison
selSamples <- unique(relativePerformanceTable[relativePerformanceTable$titrationType%in%selSamples & relativePerformanceTable$region%in%selRegions & relativePerformanceTable$dataset%in%selDatasets & relativePerformanceTable$expDifference>1, "titrationType"])

perfTableRegs <- unique(relativePerformanceTable$region)
perfTableCmps <- unique(relativePerformanceTable$titrationType)

# compare each dataset against each other datasets compiling a 3x3 contigency/concordance table:
allTabs <- list()
relativeAssayVsAssayData <- data.frame()
relativeAssayVsAssayDataConcord <- matrix(NA,nrow=length(selDatasets),ncol=length(selDatasets),dimnames=list(selDatasets,selDatasets))
for(ds1 in selDatasets) {
	allTabsTmp <- list()

	# get the data for the first dataset to be compared (focusing only on the relative change -1, 0, or 1):
	d1 <- relativePerformanceTable[relativePerformanceTable$dataset==ds1 & relativePerformanceTable$titrationType%in%selSamples & relativePerformanceTable$region%in%selRegions,]
	dComplete1 <- matrix(NA,ncol=length(perfTableCmps),nrow=length(perfTableRegs),dimnames=list(perfTableRegs,perfTableCmps))
	for(i in 1:nrow(d1)) dComplete1[d1$region[i],d1$titrationType[i]] <- as.numeric(d1$relativeChange[i])
	dComplete1 <- factor(as.vector(dComplete1),levels=as.character(c(-1,0,1)))

	for(ds2 in selDatasets) {
		# get the data for the second dataset to be compared (focusing only on the relative change -1, 0, or 1):
		d2 <- relativePerformanceTable[relativePerformanceTable$dataset==ds2 & relativePerformanceTable$titrationType%in%selSamples & relativePerformanceTable$region%in%selRegions,]
		dComplete2 <- matrix(NA,ncol=length(perfTableCmps),nrow=length(perfTableRegs),dimnames=list(perfTableRegs,perfTableCmps))
		for(i in 1:nrow(d2)) dComplete2[d2$region[i],d2$titrationType[i]] <- as.numeric(d2$relativeChange[i])
		dComplete2 <- factor(as.vector(dComplete2),levels=as.character(c(-1,0,1)))
	
		# tabularize the values in both datasets: how often are both 1, both 0, only the first 0, etc.?
		tab <- table(dComplete1,dComplete2)
		tab <- tab / sum(tab)
		allTabsTmp[[ds2]] <- tab
		melted <- melt(tab)
		
		# calculate the "concordance" as the proportion of the values on the diagonal divided by the sum of all values:
		relativeAssayVsAssayDataConcord[ds1,ds2] <- sum(diag(tab),na.rm=T)/sum(tab,na.rm=T)
		relativeAssayVsAssayData <- rbind(relativeAssayVsAssayData, cbind(ds1=ds1,ds2=ds2,melted))
	}

	allTabs[[ds1]] <- allTabsTmp
}

# plot the results in a symmetric heatmap (Figure 3B):

message("\t* agreement of relative assays with absolute assays (Figure 3.b)")
ggData <- relativeAssayVsAssayDataConcord
rownames(ggData) <- datasetTable[rownames(ggData),"prettyLabel"]
colnames(ggData) <- datasetTable[colnames(ggData),"prettyLabel"]
clust <- hclust(dist(ggData),"average")
svgPlot("clust_chisq_relative", 5, 5)
par(mar=c(4,2,1,8))
dend <-as.dendrogram(clust)
plot(dend,horiz=T)
dev.off()
cOrder <- rank(clust$order)
svgPlot("clust_chisq_heatmap_relative", nrow(ggData)/2, nrow(ggData)/2)
pheatmap(ggData[cOrder,cOrder],col=colorRampPalette(rev(brewer.pal(5,"Purples")),bias=0.6,space="Lab")(20),breaks=seq(0,1,0.05),border_color="white",cellwidth=22,cellheight=22,cluster_cols=FALSE,cluster_rows=FALSE,display_numbers=TRUE,fontsize=13)
dev.off()


# also plot a matrix of concordance tables to give a more detailed view of how the relative assays relate to each other (Figure 3A):
message("\t* agreement of relative assays with each other (Figure 3.a)")
cols <- c("snow","maroon")
ggData <- relativeAssayVsAssayData[relativeAssayVsAssayData$ds1 %in% datasetsByType$relative & relativeAssayVsAssayData$ds2 %in% datasetsByType$relative, ]
ggData$dComplete1 <- factor(ggData$dComplete1,levels=c(-1,0,1),labels=c("-","=","+"))
ggData$dComplete2 <- factor(ggData$dComplete2,levels=c(-1,0,1),labels=c("-","=","+"))
ggData$Pretty_num <- sprintf("%.2f", ggData$value)
d <- ggplot(ggData, aes(dComplete1, dComplete2)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low=cols[1], high=cols[2], na.value="white", guide=guide_colorbar(title="Count",bin=4)) + geom_text(aes(label=Pretty_num),family="Arial",colour="black",size=2.6) + xlab(NULL) + ylab(NULL) + coord_flip() + facet_grid(ds1 ~ ds2, scales="fixed", margins=F) + theme_bw() + theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			panel.margin = unit(0, "lines"),
			axis.title= element_text(size=8,face="bold"),
			axis.text.y = element_text(hjust=1,vjust=0.5,size=8),
			axis.text.x = element_text(hjust=1,vjust=0.5,size=8),
			axis.ticks = element_blank(), strip.background = element_blank(), strip.text = element_text(size=8,face="bold")
	) 
svgPlotGG(d, "hm_rel_assay_cmp", 14, 12)





save.image(paste("methBench","afterFig3",format(Sys.time(), "%Y%m%d_%H%M"),"session.Rdata",sep="_"))
