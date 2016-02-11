# SUPPLEMENTARY FIGURES, TABLES, AND STATISTICS


# ensure reproducibility by fixing the random seed to an arbitrary offset:
load(latestPrepName)
set.seed(1234)
message("=== SUPPLEMENTARY FIGURE AND DATA ===")




##### HEATMAP-LIKE OVERVIEW OF AVAILABLE DATASETS (Fig. S1) #####


message("Data overview table (Figure S1)...")

# calculate an overview of available assays:
tmp <- lsaData

# There where some region/sample combinations that were successfully measured by MS-HRM, 
# but which did not result in valid numeric values in the way we interpreted them. In order
# not to disadvantage MS_HRM in the suppementary figure, we fill in dummy values for these 
# measurements (but only for the figure):
tmp <- rbind(tmp,data.frame(
	"regionName"=c("region_02","region_02","region_02","region_02","region_02","region_02","region_02","region_02","region_02","region_02","region_06","region_06","region_14","region_14","region_15","region_15","region_15","region_15"),
	"sampleName"=c("KG1a_AzaC","KG1a_Control","Titration1_1_100","Titration1_2_75","Titration1_3_50","Titration1_4_10","Titration2_1_100","Titration2_2_10","Xenograft1_FFPE","Xenograft1_Frozen","CRC_1_Control","CRC_2_Control","Titration1_3_50","Titration1_5_1","Titration2_3_1","Titration2_4_01","Titration2_5_001","Titration2_6_0"),
	"combinedName"=NA,
	"methValue"=-1,
	"labelBySample"=NA,
	"labelByLocus"=NA,
	"labelCombined"=NA,
	"datasetName"="MS_HRM",
	"min"=NA, "max"=NA, "median"=NA, "mean"=NA, "lower"=NA, "upper"=NA, "deviationMedian"=NA, "deviationCorridor"=NA
))

availAssays <- aggregate(methValue~datasetName+sampleName,tmp,function(x){sum(!is.na(x))})
colnames(availAssays) <- c("datasetName","sampleName","regionCount")
availAssays$sampleName <- factor(availAssays$sampleName,levels=(sampleNames))
availAssays$datasetName <- factor(datasetTable[availAssays$datasetName,"prettyLabel"],levels=rev(datasetTable[datasetNamesInOrder,"prettyLabel"]))
availAssays$prettyNum <- sprintf("%d", availAssays$regionCount)
availAssays <- aggregate(regionName~datasetName+sampleName,tmp,function(x){ c("core"=length(intersect(regionNames[1:15],as.character(x))),"additional"=length(intersect(regionNames[16:length(regionNames)],as.character(x)))) })
availAssays <- cbind(availAssays[,1:2],availAssays$regionName,"total"=rowSums(availAssays$regionName))
availAssays$pretty <- sprintf("%d", availAssays$total)
nAdditional <- length(regionNames)-15

# assign colors:
plotColorCategories <- c(
	"All core regions, no additional regions",
	"Core regions missing, no additional regions",
	"All core regions, some additional regions",
	"Core regions missing, some additional regions",
	"All core regions, many additional regions",
	"Core regions missing, many additional regions"
)
availAssays$color <- factor((availAssays$core==15) + floor((availAssays$additional+(nAdditional/2-1)) / (nAdditional/2))*10, levels=rev(c(1,0,11,10,21,20)),labels=rev(plotColorCategories))
availAssays$sampleName <- factor(availAssays$sampleName,levels=(sampleNames))
availAssays$datasetName <- factor(datasetTable[availAssays$datasetName,"prettyLabel"],levels=rev(datasetTable[datasetNamesInOrder,"prettyLabel"]))

# generate the plot:
fillCols <- c("#777777","#cccccc","#1f78b4","#a6cee3","#6a3d9a","#cab2d6")
names(fillCols) <- plotColorCategories
textCols <- c("white","black","white","black","white","black")
names(textCols) <- plotColorCategories
d <- ggplot(availAssays, aes(sampleName, datasetName)) + geom_tile(aes(fill = color), colour = "white") + scale_fill_manual(values=fillCols) + geom_text(aes(label=pretty, colour = color),family="Arial",size=1.4) + scale_colour_manual(values=textCols) + xlab(NULL) + ylab(NULL) + defaultPlotTheme(flipX=TRUE)
svgPlotGG(d, "hm_regions_per_assay_sample", 26, 12)



### MDS PLOTS (Fig. Sx) ###

# generate a panel of MDS plots per dataset as a quick overview of the data:
message("Multi-dimensional scaling plots (Figure Sx)...")

svgPlot("mds_plots", 12, 10, pointsize=13)
par(mar=c(3,3,4,1))
par(mfrow=c(4,6))
for(curDataset in datasetNames) {

	# get current data:
	curData <- sapply(sampleNames,function(y) { allMeth[allMeth$sample==y,curDataset] })
	curData <- curData[rowSums(!is.na(curData))>0,]
	curData <- curData[,colSums(!is.na(curData))>0]

	# calculate the Euclidean distance matrix:
	d <- dist(t(curData))

	# perform multi-dimensional scaling:
	loc <- cmdscale(d)
	x <- loc[, 1]
	y <- -loc[, 2]

	# determine the geometric mean per group of the points in the plot:
	# (for label positioning)
	xGroupMedians <- aggregate(x,list(as.factor(sampleTable[rownames(loc),"sampleType"])),median,na.rm=TRUE)
	yGroupMedians <- aggregate(y,list(as.factor(sampleTable[rownames(loc),"sampleType"])),median,na.rm=TRUE)
	rownames(yGroupMedians) <- yGroupMedians[,1]

	# plot the results:
	plot(x, y, xlab = "", ylab = "", asp = 1, axes = TRUE, main = datasetTable[curDataset,"prettyLabel"], pch=16, bty="l", col=plotColLookup$sampleName[rownames(loc)])

	# add labels at the center of each group:
	par(xpd=TRUE)
	text(xGroupMedians[,2],yGroupMedians[xGroupMedians[,1],2],xGroupMedians[,1],col=plotColLookup$sampleType[as.character(xGroupMedians[,1])])	
	par(xpd=FALSE)

}
# add a legend:
plot.new()
par(xpd=NA)
legend("top",names(plotColLookup$sampleType),bty="n",border=NA,fill=plotColLookup$sampleType,cex=1.2)
dev.off()




### FULL PANEL OF BEESWARM PLOTS (Fig. S4) ###

message("Beeswarm plots panel for all assays (Figure S4.a)...")

selRegions <- coreRegionNames
selSamples <- setdiff(sampleNames,sampleNames[grep("Titration",sampleNames)])
for(curRegion in selRegions) {
	plot3WayBeeswarm(lsaData[ lsaData$regionName==curRegion  & lsaData$sampleName%in%selSamples & lsaData$datasetName%in%datasetsByType$absolute,], extension=curRegion, meth="swarm", bee.cex=1.1, bee.spacing=1, bee.main=curRegion)
}

# determine which datasets define each consensus corridor:

message("Contribution to consensus heatmap (Figure S4.b)...")

#selSamples <- sampleNames
inConsensus <- aggregate(assayGroup~datasetName+regionName+sampleName,consensusFreq[consensusFreq$sampleName%in%selSamples & consensusFreq$regionName%in%selRegions,],function(x){ 1 })
inConsensus$datasetName <- factor(datasetTable[inConsensus$datasetName,"prettyLabel"],levels=rev(datasetTable[datasetNamesInOrder,"prettyLabel"]))
colnames(inConsensus) <- c("datasetName","regionName","sampleName","color")
inConsensus$color <- plotColLookup$prettyLabel[as.character(inConsensus$datasetName)]

tmpMeths <- is.na(allMeth[,datasetsByType$absolute]) & allMeth$region%in%selRegions & allMeth$sample%in%selSamples
for(curDataset in datasetsByType$absolute) {
	naDataPoints <- allMeth[tmpMeths[,curDataset],c("sample","region")]
	if(nrow(naDataPoints)>0) {
		inConsensus <- rbind(inConsensus,data.frame(datasetName=datasetTable[curDataset,"prettyLabel"], regionName=naDataPoints$region, sampleName=naDataPoints$sample, color="gray"))
	}
}
d <- ggplot(inConsensus, aes(paste(sampleName,regionName), datasetName)) + geom_tile(aes(fill = color), colour = "white") + xlab("Region/sample combination") + ylab(NULL) + scale_fill_identity(guide = FALSE) + theme_bw() + theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_rect(colour="black"),
			axis.title= element_text(size=8,face="bold"),
			axis.text.y = element_text(hjust=1,vjust=0.5,size=8),
			axis.text.x = element_blank(),
			axis.ticks = element_blank()
	) 
svgPlotGG(d,"contribution_to_consensus", 14, 6)
# the plotted stripes in the heatmap are often dropped during conversion (on journal submission, etc.), so we also plot them in a rendered version and swap those into the respective part of the figure panel:
ggsave(filename=paste("results_analysis/plots/contribution_to_consensus.png",sep=""), plot=d, scale=1, width=14*4, height=6*4, units="cm")

# summarize the contributing dataset per assay type:
tmp <- merge(
	aggregate(regionName ~ assayGroup, consensusFreq[consensusFreq$sampleName%in%selSamples & consensusFreq$regionName%in%selRegions,], length),
	aggregate(methValue ~ assayGroup, cbind("methValue"=lsaData[lsaData$datasetName%in%datasetsByType$absolute & lsaData$sampleName%in%selSamples & lsaData$regionName%in%selRegions,"methValue"],"assayGroup"=datasetTable[lsaData[lsaData$datasetName%in%datasetsByType$absolute  & lsaData$sampleName%in%selSamples & lsaData$regionName%in%selRegions,"datasetName"],"assayGroup"]),function(x) { sum(!is.na(x)) }),
	by="assayGroup"
)
colnames(tmp) <- c("Assay type","Constituent of consensus corridor","Not contributing to corridor")
tmp[,"Not contributing to corridor"] <- tmp[,"Not contributing to corridor"] - tmp[,"Constituent of consensus corridor"]
tmp <- melt(tmp)
tmp$percent <- NA
tmp$height <- NA
tmp$height[tmp$variable=="Constituent of consensus corridor"] <- tmp[tmp$variable=="Constituent of consensus corridor","value"]+tmp[tmp$variable=="Not contributing to corridor","value"]
tmp$percent[tmp$variable=="Constituent of consensus corridor"] <- sprintf("%.1f %%", tmp[tmp$variable=="Constituent of consensus corridor","value"] / tmp$height[tmp$variable=="Constituent of consensus corridor"] * 100.0) 
tmp$height <- tmp$height+100
tmp$color <- plotColLookup$assayGroup[as.character(tmp[,"Assay type"])]
tmp$color[tmp$variable=="Not contributing to corridor"] <- makeTransparent(tmp$color[tmp$variable=="Not contributing to corridor"],0.25) 

# and visualize the summary as a barplot:
message("Contribution to consensus bar plots (Figure S4.c)...")
d <- ggplot(tmp, aes(x=`Assay type`,y=value,fill=color)) + geom_bar(width=0.5,colour="#333333",size=0.25,stat="identity") + geom_hline(yintercept=0) + ylim(0,1.2*max(tmp$height,na.rm=TRUE))+ xlab(NULL) + geom_text(aes(y=height,label=percent),angle=90,size=2,hjust=0) + ylab("Number of measurements") + scale_fill_identity() + defaultPlotTheme(flipX=TRUE,fontSize=8)  
svgPlotGG(d, "contribution_to_consensus_summary", 4, 6)


### FRESH-FROZEN VS FFPE ###


message("Fresh-frozen vs FFPE (Figures Sxxx)...")

selSamples <- sampleNamesByType$FrozenFFPE
selRegions <- regionNames
selDatasets <- datasetsByType$absolute

# boxplots for absolute deviation to corridor:
message("\t* absolute deviation")
ggData <- melt(lsaData[lsaData$datasetName%in%selDatasets & lsaData$regionName%in%selRegions & lsaData$sampleName%in%selSamples,][,c("datasetName","deviationCorridor","sampleName")],id.vars=c("datasetName","sampleName"))
ggData$sampleType <- ifelse(grepl("Frozen",ggData$sampleName),"Frozen","FFPE")
ggData <- cbind(ggData,"col"=makeTransparent(plotColLookup$datasetName[ggData$datasetName],0.75))
ggData$datasetNameAndSampleType <- factor(paste(datasetTable[ggData$datasetName,"prettyLabel"],ggData$sampleType),levels=c(sapply(datasetTable[datasetNamesInOrder,"prettyLabel"], paste, c("Frozen","FFPE"))))
ggData$variable <- "Relative to consensus corridor"
d <- ggplot(ggData, aes(x=datasetNameAndSampleType, y=abs(value), fill=col)) + geom_boxplot(notch=FALSE,outlier.size=1.5,type="count") +ylim(0,100) + scale_fill_identity(guide = FALSE) + xlab(NULL) + ylab("Absolute deviation") + defaultPlotTheme(flipX=TRUE)
svgPlotGG(d,"performance_frozenffpe", 18.3, 12, units="cm")

# boxplots for "bias" (i.e. trend to over- or under-estimate value):
message("\t* directional deviation")
ggData2 <- melt(lsaData[lsaData$datasetName%in%selDatasets & lsaData$regionName%in%selRegions & lsaData$sampleName%in%selSamples,][,c("datasetName","deviationCorridor","sampleName")],id.vars=c("datasetName","sampleName"))
ggData2$sampleType <- ifelse(grepl("Frozen",ggData$sampleName),"Frozen","FFPE")
ggData2 <- aggregate(value ~ datasetName + variable + sampleType, ggData2, mean, simplify=T)
ggData2 <- cbind(ggData2,"col"=makeTransparent(plotColLookup$datasetName[ggData2$datasetName],0.75))
ggData2$datasetNameAndSampleType <- factor(paste(datasetTable[ggData2$datasetName,"prettyLabel"],ggData2$sampleType),levels=c(sapply(datasetTable[datasetNamesInOrder,"prettyLabel"], paste, c("Frozen","FFPE"))))
ggData$variable <- "Relative to consensus corridor"
d <- ggplot(ggData2, aes(x=datasetNameAndSampleType,y=value,fill=value)) + geom_bar(stat="identity",colour="#333333") + ylim(-5,5) + xlab(NULL) + ylab("Bias")  + scale_fill_gradient2(low=brewer.pal(3,"BrBG")[3],high=brewer.pal(3,"BrBG")[1],guide=F,space="Lab") + geom_hline(aes(yintercept=0),colour="#333333") + defaultPlotTheme(flipX=TRUE)  
svgPlotGG(d,"performance_frozenffpe_bias_quant", 18.3, 9, units="cm")


# version 2: always compare to "frozen" corridor (even for FFPE):

lsaDataMod <- lsaData[lsaData$datasetName%in%setdiff(selDatasets,"Infinium") & lsaData$regionName%in%selRegions & lsaData$sampleName%in%selSamples,]

for(curRegion in selRegions) {
	for(curSample in selSamples) {
		if(grepl("FFPE",curSample)) {
			replWithSample <- gsub("FFPE","Frozen",curSample)
			for(curDataset in selDatasets) {
				repl <- lsaDataMod[lsaDataMod$sampleName==replWithSample & lsaDataMod$datasetName==curDataset & lsaDataMod$regionName==curRegion,c("upper","lower","median")]
				if(nrow(repl)>0) {
					lsaDataMod[lsaDataMod$sampleName==curSample & lsaDataMod$datasetName==curDataset & lsaDataMod$regionName==curRegion,c("upper","lower","median")] <- repl
				}
				else {
					lsaDataMod[lsaDataMod$sampleName==curSample & lsaDataMod$datasetName==curDataset & lsaDataMod$regionName==curRegion,c("upper","lower","median")] <- NA
				}
			}
		}
	}
}
lsaDataMod$deviationMedian <- lsaDataMod$methValue - lsaDataMod$median
lsaDataMod$deviationCorridor <- ifelse(lsaDataMod$methValue>=lsaDataMod$lower & lsaDataMod$methValue<=lsaDataMod$upper,0,ifelse(lsaDataMod$methValue<lsaDataMod$lower, lsaDataMod$methValue-lsaDataMod$lower, lsaDataMod$methValue-lsaDataMod$upper))

# boxplots for absolute deviation to corridor:
message("\t* absolute deviation")
ggData <- melt(lsaDataMod[,c("datasetName","deviationCorridor","sampleName")],id.vars=c("datasetName","sampleName"))
ggData$sampleType <- ifelse(grepl("Frozen",ggData$sampleName),"Frozen","FFPE")
ggData <- cbind(ggData,"col"=makeTransparent(plotColLookup$datasetName[ggData$datasetName],0.75))
ggData$datasetNameAndSampleType <- factor(paste(datasetTable[ggData$datasetName,"prettyLabel"],ggData$sampleType),levels=c(sapply(datasetTable[datasetNamesInOrder,"prettyLabel"], paste, c("Frozen","FFPE"))))
ggData$variable <- "Relative to consensus corridor"
d <- ggplot(ggData, aes(x=datasetNameAndSampleType, y=abs(value), fill=col)) + geom_boxplot(notch=FALSE,outlier.size=1.5,type="count") +ylim(0,100) + scale_fill_identity(guide = FALSE) + xlab(NULL) + ylab("Absolute deviation") + defaultPlotTheme(flipX=TRUE)
svgPlotGG(d,"performance_frozenffpe_v2", 18.3, 12, units="cm")

# boxplots for "bias" (i.e. trend to over- or under-estimate value):
message("\t* directional deviation")
ggData2 <- melt(lsaDataMod[,c("datasetName","deviationCorridor","sampleName")],id.vars=c("datasetName","sampleName"))
ggData2$sampleType <- ifelse(grepl("Frozen",ggData$sampleName),"Frozen","FFPE")
ggData2 <- aggregate(value ~ datasetName + variable + sampleType, ggData2, mean, simplify=T)
ggData2 <- cbind(ggData2,"col"=makeTransparent(plotColLookup$datasetName[ggData2$datasetName],0.75))
ggData2$datasetNameAndSampleType <- factor(paste(datasetTable[ggData2$datasetName,"prettyLabel"],ggData2$sampleType),levels=c(sapply(datasetTable[datasetNamesInOrder,"prettyLabel"], paste, c("Frozen","FFPE"))))
ggData$variable <- "Relative to consensus corridor"
d <- ggplot(ggData2, aes(x=datasetNameAndSampleType,y=value,fill=value)) + geom_bar(stat="identity",colour="#333333") + ylim(-12,5) + xlab(NULL) + ylab("Bias")  + scale_fill_gradient2(low=brewer.pal(3,"BrBG")[3],high=brewer.pal(3,"BrBG")[1],guide=F,space="Lab") + geom_hline(aes(yintercept=0),colour="#333333") + defaultPlotTheme(flipX=TRUE)  
svgPlotGG(d,"performance_frozenffpe_bias_quant_v2", 18.3, 9, units="cm")



### FULL SCATTER PLOT PANELS OF ABSOLUTE AND RELATIVE ASSAYS (Supplementary Fig. 2 and 8) ###

selSamples <- sampleNames
message("Scatter plot panels (Supplementary Figures 2 and 8)...")

scatterComparisons <- list("absVsAbs"=c("absolute","absolute"),"relVsAll"=c("relative","absoluteAndRelative"))
for(scatterCmp in scatterComparisons) {
	# select data:
	selDatasets1 <- datasetsByType[[scatterCmp[1]]]
	selDatasets2 <- datasetsByType[[scatterCmp[2]]]
	selIndices1 <- assayVsAssayData$assay1%in%selDatasets1
	selIndices2 <- assayVsAssayData$assay2%in%selDatasets2
	selDatasets <- union(selDatasets1,selDatasets2)
	ggData <- assayVsAssayData[selIndices1 & selIndices2,]
	ggData$assay1 <- factor(datasetTable[as.character(ggData$assay1),"prettyLabel"],levels=datasetTable[selDatasets,"prettyLabel"])
	ggData$assay2 <- factor(datasetTable[as.character(ggData$assay2),"prettyLabel"],levels=datasetTable[selDatasets,"prettyLabel"])

	# get correlation table (pre-computed!): 
	corTablePretty <- corTable
	rownames(corTablePretty) <- datasetTable[rownames(corTablePretty),"prettyLabel"]
	colnames(corTablePretty) <- datasetTable[colnames(corTablePretty),"prettyLabel"]

	# augment data table by adding the correlation coefficient numbers as "pseudo-rows"
	# to the original table of DNA methylation measurements
	a1 <- unique(ggData[as.numeric(ggData$assay1) <= as.numeric(ggData$assay2),]$assay1)
	a2 <- unique(ggData[as.numeric(ggData$assay1) <= as.numeric(ggData$assay2),]$assay2)
	ggData <- ggData[(as.numeric(ggData$assay1) < as.numeric(ggData$assay2)),] # upper half only
	ggData$sample <- as.character(ggData$sample)
	ggData$txt <- NA
	ggData$txtX <- 50
	ggData$txtY <- 50
	ggData$txtCol <- 50
	for(s1 in a1) {
		for(s2 in a2) {
			# doesn't work due to bug in ggplots --> set text to transparent instead:
			#if(which(s1==levels(ggData$assay1))<which(s2==levels(ggData$assay2))) next;
			corVal <- corTablePretty[s1,s2]
			txt <- ifelse(s1==s2,s1,sprintf("r = %.2f", corVal)) 
			txtX <- 50
			txtY <- 50
			txtCol <- "black"
			if(which(s1==levels(ggData$assay1))<which(s2==levels(ggData$assay2))) {
				txtCol <- rgb(0,0,0,alpha=0)
				#if(scatterCmp[1]==scatterCmp[2]) {
					# make redundant correlation values transparent (otherwise they'll overlap with points!) 
				#	txtCol <- rgb(0,0,0,alpha=0)
				#}
				#else {
				#	txtX <- 80
				#	txtY <- 3
				#}
			}
			ggData <- rbind(ggData,list("assay1"=s1,"assay2"=s2,"region"="none","sample"="none","meth1"=NA,"meth2"=NA,"txt"=txt,"txtX"=txtX,"txtY"=txtY,"txtCol"=txtCol))
		}
	}

	# generate a scatterplot panel:
	d <- ggplot(data=ggData,aes(x=meth1,y=meth2)) + geom_point(shape=16,alpha=0.3)  + geom_smooth(method=lm,fill='steelblue') + geom_text(aes(label=txt,x=txtX,y=txtY,colour=txtCol), size=5, stat = "identity", position = "identity") + facet_grid(assay2 ~ assay1 ,margins=FALSE,drop=T) +scale_colour_identity() + xlab(DNA_METH_LABEL) + ylab(DNA_METH_LABEL) + xlim(0,100) + ylim(0,100)  + theme_bw() + theme(	axis.line = element_line(), panel.grid.major = element_blank(),	panel.grid.minor = element_blank(), axis.title= element_text(size=16,face="bold"), axis.text.y = element_text(hjust=1,vjust = 0.5,size = 8),	axis.text.x = element_text(size = 8), strip.background = element_blank(), strip.text = element_text(size=14,face="bold")) 
	svgPlotGG(d, paste0("pairwise_scatterplots_ggplot_matrix_",scatterCmp[1],"-vs-",scatterCmp[2]), 1.65*length(selDatasets1), 1.6*length(selDatasets2), units="in") 

	# plot correlation heatmap
	svgPlot(paste0("pairwise_scatterplots_ggplot_matrix_cor_heatmap_",scatterCmp[1],"-vs-",scatterCmp[2]), 12, 12)
	pheatmap(corTablePretty,col=colorRampPalette(rev(brewer.pal(5,"Oranges")),bias=0.2,space="Lab")(20),breaks=seq(0,1,0.05),number_color="black",border_color="white",cellwidth=22,cellheight=22,cluster_cols=FALSE,cluster_rows=FALSE,display_numbers=TRUE,fontsize=13,fontcolor="black")
	dev.off()

}



##### LINEAR MODELS FOR CAUSES OF DEVIATION (Supplementary Figure 4) #####

message("Causes of deviation (Supplementary Figure 4)...")

message("\t* build models")

# choose data to work with:
selDatasets <- datasetsByType$absolute
selSamples <- sampleNames
dependentVar <- "deviationCorridor" # this is variable that we'd like to find explanations for

# annotate lsaData with additional information (from the region annotation table):
annotatedLsaData <- merge(lsaData,regionTable,by="regionName")
annotatedLsaData <- annotatedLsaData[order(annotatedLsaData$regionName),]
annotatedLsaData <- annotatedLsaData[order(annotatedLsaData$sampleName),]
annotatedLsaData <- annotatedLsaData[annotatedLsaData$sampleName %in% selSamples & annotatedLsaData$datasetName %in% selDatasets, ]
annotatedLsaData[,dependentVar] <- abs(annotatedLsaData[,dependentVar]) # absolute deviation!

# examine the following columns as potential causes of deviations:
selColsFactor <- c("datasetName") 
selColsNumeric <- c("median","gCContent","obsExpRatio","repeatContent") 
selColsBinary <- character(0)

# convert factors into binary indicator variables:
for (curCol in selColsFactor) {
  a <- model.matrix(as.formula(paste("~",curCol,"-1",sep="")), data=annotatedLsaData)
  colnames(a) <- gsub(curCol,"",colnames(a))
  selColsBinary <- c(selColsBinary,colnames(a))
  annotatedLsaData <- cbind(annotatedLsaData,a)
}
selColsAll <- c(selColsNumeric,selColsBinary)

# standardize numeric values into an interval from zero to one:
for (curCol in selColsNumeric) {
  curVals <- annotatedLsaData[,curCol]
  minVal <- unname(quantile(curVals,0.005,na.rm=T))
  maxVal <- unname(quantile(curVals,0.995,na.rm=T))
  annotatedLsaData[,curCol] = (annotatedLsaData[,curCol]-minVal)/(maxVal-minVal)
}   

# define the linear model to be examined, i.e. deviation from consensus corridor
# as a function of dataset and (potential) explanatory region characteristics:
curModel <- as.formula(paste(
	dependentVar, "~"
	,"(",paste(selDatasets,collapse=" + "),")"
	,"*"
	,"(",paste(selColsNumeric,collapse=" + "),")"
	,"-"
	,paste(selColsNumeric,collapse=" - ")
))

# fit the model and calculate p-values:
modelData <- annotatedLsaData[annotatedLsaData$regionName %in% regionNames,c("regionName","sampleName","datasetName",dependentVar,selColsAll)]
modelData$datasetName <- as.factor(modelData$datasetName)
fit <- lm(formula=curModel,data=modelData)
sfit <- summary(fit)  

# for doing some variable re-labeling:
recodeVarType <- function(x) {
	factor(recode(as.character(x),"'methValue'='DNA methylation (%)';'median'='Consensus corridor (median)';'gCContent'='GC content';'obsExpRatio'='CpG observed/expected ratio';'repeatContent'='Repetitive DNA content'"),levels=c("DNA methylation (%)","Consensus corridor (median)","GC content","CpG observed/expected ratio","Repetitive DNA content"))
}

# summarize modeling results in a table:
lmResults <- numeric()
for(curDataset in unique(modelData$datasetName)) {
	for(curVar in selColsNumeric) {		
		curName <- paste(curDataset,":",curVar,sep="")

		curVec <- list("dataset"=curDataset, "variable"=curVar, "dir"=NA, "p"=NA, "stepDir"=NA, "stepP"=NA)		

		if(curName %in% rownames(sfit$coefficients)) {
			tmp <- as.numeric(sfit$coefficients[curName,])
			curVec$dir <- sign(tmp[3])
			curVec$p <- tmp[4]
		}
		
		lmResults <- rbind(lmResults,curVec)
	}
}
lmResults <- as.data.frame(lmResults)
rownames(lmResults) <- NULL
lmResults$dir <- as.numeric(as.character(lmResults$dir))
lmResults$p <- as.numeric(as.character(lmResults$p))
lmResults$adjP <- p.adjust(lmResults$p,method="BH") # adjust all p-values for multiple hypothesis testing!
lmResults$mLog10P <- -log10(lmResults$adjP) * lmResults$dir # this is a transformation of the p-value so that smaller p-values get larger numbers (larger bars in the plot later) and so the negative values represent an inverse relationship and positive values a direct relationship between factor and effect
lmResults$color <- plotColLookup$datasetName[as.character(lmResults$dataset)] 
lmResults$variable <- recodeVarType(lmResults$variable)
lmResults$dataset <- factor(datasetTable[as.character(lmResults$dataset),"prettyLabel"],levels=(datasetTable[datasetNamesInOrder,"prettyLabel"]))

ggData <- lmResults

# the plot contains a few very extreme p-values, which need to be cropped so that the other p-values can be
# better visualized. We therefore trim mLog10P values above or below the 95% quantile (these values will be 
# annotated in the plot later on):
maxP <- ceiling(quantile(abs(ggData$mLog10P),0.95,na.rm=T)+1)
ggData$mLog10P[ggData$mLog10P > maxP] <- maxP
ggData$mLog10P[ggData$mLog10P < -maxP] <- -maxP

message("\t* make plots")

# generate a panel of barplots:
d <- ggplot(ggData, aes(x=dataset,y=mLog10P,fill=color)) + geom_bar(width=0.5,colour="#333333",size=0.25,stat="identity")  + scale_fill_identity() + geom_hline(yintercept=0) + geom_hline(yintercept=+log10(0.05),colour="gray",linetype=2) + geom_hline(yintercept=-log10(0.05),colour="gray",linetype=2) + ylim(c(-maxP,maxP)) + facet_wrap( ~ variable, scales="free_y", ncol=1, drop=T) + defaultTheme(flipX=TRUE) + xlab(NULL)
suppressWarnings(ggsave(filename=paste("results_analysis/plots/influences.svg",sep=""), plot=d, scale=1, width=11, height=16, units="cm"))

# generate illustrative scatter plots for those assays that are most strongly influenced by each region characteristic:
ggScatterData <- melt(modelData[,c("datasetName",dependentVar,selColsNumeric)],id.vars=c("datasetName",dependentVar))
ggScatterData$datasetName <- factor(datasetTable[as.character(ggScatterData$datasetName),"prettyLabel"],levels=rev(datasetTable[datasetNamesInOrder,"prettyLabel"]))
for(curVar in selColsNumeric) {	
	tmp <- ggData[ggData$variable==recodeVarType(curVar),]
	tmp <- tmp[tmp$adjP==min(tmp$adjP),]
	tmp2 <- ggScatterData[ggScatterData$datasetName==tmp$dataset & ggScatterData$variable==curVar,]
	maxY <- 100

	svg(filename=paste("results_analysis/plots/scatter_influences_",curVar,".svg",sep=""),width=2.4, height=2.6)
	par(cex=1)
	par(mar=c(4,4,2,1))
	plot(tmp2$value, tmp2$deviationCorridor, pch=16, xlab=recodeVarType(curVar), ylab="Abs. dev. from corridor", las=1,col=rgb(0,0,0,0.3), ylim=c(0,maxY), xlim=c(0,1), main=tmp$dataset)
	abline(lm(deviationCorridor ~ value, tmp2), col="hotpink")
	par(xpd=NA)
	text(0.5,y=1.1 * maxY,sprintf("r = %.1f, adj.P = %s", cor(tmp2$deviationCorridor, tmp2$value, use="pairwise"), format(tmp$adjP, scientific = TRUE, digits = 2)))
	dev.off()
}




##### LINEAR MODELS FOR INFLUENCE OF DNA AMOUNT ON DEVIATION (Supplementary Figure 9) #####

message("Influence of DNA amount on deviations (Supplementary Figure 9)...")

message("\t* build models")

# choose data to work with:
selDatasets <- datasetsByType$absoluteAndRelative
selSamples <- sampleNames
dependentVar <- "deviationCorridor" #"deviationCorridor" # this is variable that we'd like to find explanations for

# load sample quantities data:
sampleQuantities <- read.table("data/sample_quantities.txt",sep="\t",header=TRUE)
sampleQuantities$dataset <- gsub(" ","_",sampleQuantities$dataset)
sampleQuantities$sample <- sampleNameLookup[sanitizeSampleNames(sampleQuantities$sample)]

for(multiSample in unique(grep(",",sampleQuantities$dataset,value=TRUE))) {
	curDatasets <- strsplit(multiSample,",")[[1]]
	sampleQuantities[sampleQuantities$dataset==multiSample,"dataset"] <- curDatasets[1]
	for(ds in curDatasets[-1]) {
		tmp <- sampleQuantities[sampleQuantities$dataset==curDatasets[1],]
		tmp$dataset <- ds
		sampleQuantities <- rbind(sampleQuantities,tmp)
	}
}
rownames(sampleQuantities) <- paste(sampleQuantities$dataset,sampleQuantities$sample,sep="_")

sampleQuantities$measurementTech <- as.factor(tolower(sapply(strsplit(as.character(sampleQuantities$measurementTaken)," "),function(x) x[1])))
sampleQuantities$measurementTech <- recode(sampleQuantities$measurementTech,"'picogreen'='Quant-iT';'quant-it'='Quant-iT';'qubit'='Qubit';'nanodrop'='NanoDrop';'tapestation'='TapeStation'")
measurementTechCols <- structure(brewer.pal(length(levels(sampleQuantities$measurementTech)),"Dark2"),names=levels(sampleQuantities$measurementTech))

# correlation heatmap of similarity between concentration measurements:

sampleConcTable <- dcast(sampleQuantities, sample~dataset, value.var="measuredConcentration")
rownames(sampleConcTable) <- sampleConcTable[,1]
sampleTechs <- aggregate(measurementTech~dataset,sampleQuantities,unique) 
sampleTechs <- structure(sampleTechs[,2],names=as.character(datasetTable[sampleTechs[,1],"prettyLabel"]))
sampleConcTable <- sampleConcTable[,-1]
sampleConcTable <- sampleConcTable[,colSums(!is.na(sampleConcTable))>0]
colnames(sampleConcTable) <- datasetTable[colnames(sampleConcTable),"prettyLabel"]
sampleConcTableUniq <- sampleConcTable
toMerge <- list(
	"qMSP (standard/preamp)"=c("qMSP (standard)","qMSP (preamp)"),
	"Pyroseq 1 (+replicate)"=c("Pyroseq 1","Pyroseq 1 (replicate)"),
	"Immunoquant, Pyroseq AluYb8/D4Z4/LINE1/NBL2"=c("Immunoquant","Pyroseq AluYb8","Pyroseq D4Z4","Pyroseq LINE1","Pyroseq NBL2"),
	"AmpliconBS 1, EpiTyper 1"=c("AmpliconBS 1","EpiTyper 1")
)
for(multiSample in names(toMerge)) {
	curDatasets <- toMerge[[multiSample]]
	curDatasets <- intersect(curDatasets,colnames(sampleConcTableUniq))
	if(length(curDatasets)>1) {
		sampleConcTableUniq <- sampleConcTableUniq[,setdiff(colnames(sampleConcTableUniq),curDatasets[-1])]
		colnames(sampleConcTableUniq)[colnames(sampleConcTableUniq)==curDatasets[1]] <- multiSample
		sampleTechs <- sampleTechs[setdiff(names(sampleTechs),curDatasets[-1])]
		names(sampleTechs)[names(sampleTechs)==curDatasets[1]] <- multiSample
	}
}

# plot beeswarm of input DNA concentrations:
consensusConc <- matrix(NA, ncol=2, nrow=length(sampleNames), dimnames=list(sampleNames, c("conc","relConc"))) 
plotData <- melt(data.matrix(sampleConcTableUniq))
plotData <- data.frame(plotData,sampleTechs[as.character(plotData[,2])],stringsAsFactors=FALSE)
colnames(plotData) <- c("x","dataset","y","tech")
curSamples <- sort(unique(plotData$x))
consLineWidth <- 0.25
plotData <- plotData[order(plotData$x),]
concBeePlot <- function(plotData, ylab, ylim=NA) {	
	beeswarm(y~x, data=plotData, pwcol=measurementTechCols[plotData$tech], las=3, ylab=ylab, cex=0.6, spacing=2, cex.lab=1, bty='l', cex.axis=1,method="hex",corral="none",ylim=ylim,pch=16)
}

# .. by concentration:
svgPlot("dnaconc_beeswarms", 7, 4.4)
par(mar=c(7.5,4,1,1))
concBeePlot(plotData, "Input DNA concentration [ng/ul]",ylim=c(0,440))
for(curSample in curSamples) {
	consensusConc[curSample,"conc"] <- median(plotData[plotData$x==curSample,"y"],na.rm=TRUE)
	x <- which(curSamples==curSample)
	segments(x-consLineWidth,consensusConc[curSample,"conc"],x+consLineWidth,consensusConc[curSample,"conc"],lwd=1.5,col="black")
}
par(xpd=NA)
legend("top",names(measurementTechCols),col=measurementTechCols,pch=16,horiz=TRUE,bty="n",border=NA)
dev.off()

# .. and by amount:
svgPlot("dnaconc_beeswarms_amount", 10, 4.4)
FIXED_VOL <- 25
plotDataAmount <- plotData
plotDataAmount$y <- plotDataAmount$y * FIXED_VOL
par(mar=c(7.5,4,1,1))
concBeePlot(plotDataAmount, "Amount of input DNA [ng]",ylim=c(0,440)*FIXED_VOL)
concAmounts <- matrix(NA,nrow=length(curSamples),ncol=3,dimnames=list(curSamples,c("median","mean","sd")))
for(curSample in curSamples) {
	concAmounts[curSample,] <- c(
		median(plotDataAmount[plotDataAmount$x==curSample,"y"],na.rm=TRUE),
		mean(plotDataAmount[plotDataAmount$x==curSample,"y"],na.rm=TRUE),
		sd(plotDataAmount[plotDataAmount$x==curSample,"y"],na.rm=TRUE)
	)
	x <- which(curSamples==curSample)
	segments(x-consLineWidth,concAmounts[curSample,"median"],x+consLineWidth,concAmounts[curSample,"median"],lwd=1.5,col="black")
}
par(xpd=NA)
legend("top",names(measurementTechCols),col=measurementTechCols,pch=16,horiz=TRUE,bty="n",border=NA)
dev.off()

# .. and finally the numbers:
svgPlot(paste0("dnaconc_beeswarms_amount_numbers"), 10, 2.3)
pheatmap(round(t(concAmounts)),number_color="black",border_color="white",cluster_cols=FALSE,cluster_rows=FALSE,display_numbers=TRUE,number_format="%d",fontsize=13,fontcolor="black")
dev.off()


# then model the influence of DNA amount on measurement accuracy (deviation from consensus):

modelData <- merge(lsaData[lsaData$regionName %in% coreRegionNames,],sampleQuantities,by.x=c("datasetName","sampleName"),by.y=c("dataset","sample"))
modelData <- data.frame(modelData, consensusConc[modelData$sampleName,],row.names=rownames(modelData))
modelData <- modelData[!is.na(modelData[,dependentVar]),]

# model by:
# 1. total amount used according to own measurements:
i <- is.na(modelData$measuredVolume)
modelData$measuredVolume[i] <- modelData$actualStockVolume[i]
modelData$amountUsedOwnEstimate <- modelData$measuredVolume * modelData$volumeProportionOfInput * modelData$measuredConcentration
# 2. total amount in tube according to own measurements (fixed vol 20ul):
modelData$amountUsedConsensus <- modelData$actualStockVolume * modelData$volumeProportionOfInput *  modelData$conc

# examine the following columns as potential causes of deviations:
selColsFactor <- c("datasetName","regionName") 
selColsNumeric <- c("amountUsedOwnEstimate","amountUsedConsensus") 
selColsBinary <- character(0)

# for doing some variable re-labeling:
recodeVarType <- function(x) {
	tmp <- c("Amount of input DNA by own measurement","Amount of input DNA by consensus")
	factor(recode(as.character(x),paste0("'amountUsedOwnEstimate'='",tmp[1],"';'amountUsedConsensus'='",tmp[2],"'")),levels=tmp)
}

# convert factors into binary indicator variables:
for (curCol in selColsFactor) {
  a <- model.matrix(as.formula(paste("~",curCol,"-1",sep="")), data=modelData)
  colnames(a) <- gsub(curCol,"",colnames(a))
  selColsBinary <- c(selColsBinary,colnames(a))
  modelData <- cbind(modelData,a)
}
selColsAll <- c(selColsNumeric,selColsBinary)

# perform a correlation test first:
corTestRes <- matrix(NA, ncol=length(selColsNumeric), nrow=length(unique(modelData$datasetName)), dimnames=list(unique(modelData$datasetName),selColsNumeric))
corTestDirs <- matrix(NA, ncol=length(selColsNumeric), nrow=length(unique(modelData$datasetName)), dimnames=list(unique(modelData$datasetName),selColsNumeric))
for(curDataset in unique(modelData$datasetName)) {
	curData <- modelData[modelData$datasetName==curDataset,c(dependentVar,selColsNumeric)]
	x <- abs(curData[,dependentVar])
	for(curVar in selColsNumeric) {	
		y <- curData[,curVar]
		if(sum(!is.na(y))>1 & sd(y)>0.02) { # exclude those assays without measurements or with little or variation in the estimates or real values
			tmp <- cor.test(x, y)
			corTestRes[curDataset,curVar] <- tmp$p.value
			corTestDirs[curDataset,curVar] <- sign(tmp$estimate)
		}
	}
}
corTestRes <- matrix((corTestRes),ncol=ncol(corTestRes),dimnames=dimnames(corTestRes))
corResults <- merge(melt(corTestRes),melt(corTestDirs),by=c("Var1","Var2"))
colnames(corResults) <- c("dataset","variable","p","dir")
corResults <- corResults[!is.na(corResults$p),]
corResults$adjP <- p.adjust(corResults$p,method="fdr")
corResults$mLog10P <- -log10(corResults$adjP) * corResults$dir # this is a transformation of the p-value so that smaller p-values get larger numbers (larger bars in the plot later) and so the negative values represent an inverse relationship and positive values a direct relationship between factor and effect
corResults$color <- plotColLookup$datasetName[as.character(corResults$dataset)] 
corResults$variable <- recodeVarType(corResults$variable)
corResults$dataset <- factor(datasetTable[as.character(corResults$dataset),"prettyLabel"],levels=(datasetTable[datasetNamesInOrder,"prettyLabel"]))

# define the linear model to be examined, i.e. deviation from consensus corridor
# as a function of dataset and (potential) explanatory region characteristics:
curModel <- as.formula(paste(
	dependentVar, "~"
	, paste(setdiff(selColsFactor,"datasetName"),collapse="+")
	,"+"
	,"(",paste(intersect(intersect(selDatasets,colnames(modelData)),names(which(rowSums(!is.na(corTestRes))>0))),collapse=" + "),")"	# use only the dataset with values for the measurements (determined earlier when the correlation test was calculated)
	,"*"
	,"(",paste(selColsNumeric,collapse=" + "),")"
	,"-"
	,paste(selColsNumeric,collapse=" - ")
))

# fit the model and calculate p-values:
modelData$datasetName <- as.factor(modelData$datasetName)
modelData[,dependentVar] <- abs(modelData[,dependentVar])
fit <- lm(formula=curModel,data=modelData)
sfit <- summary(fit)  

# summarize modeling results in a table:
lmResults <- numeric()
for(curDataset in unique(modelData$datasetName)) {
	for(curVar in selColsNumeric) {		
		curName <- paste(curDataset,":",curVar,sep="")

		curVec <- list("dataset"=curDataset, "variable"=curVar, "dir"=NA, "p"=NA)		

		if(curName %in% rownames(sfit$coefficients)) {
			tmp <- as.numeric(sfit$coefficients[curName,])
			curVec$dir <- sign(tmp[3])
			curVec$p <- tmp[4]
		}
		
		lmResults <- rbind(lmResults,curVec)
	}
}
lmResults <- as.data.frame(lmResults)
rownames(lmResults) <- NULL
lmResults <- lmResults[!is.na(lmResults$p) & (lmResults$variable!="Amount of input DNA by own measurement" | lmResults$dataset%in%unique(as.character(modelData[modelData$fixedByInput%in%grep("ALL",unique(modelData$fixedByInput),value=TRUE),"datasetName"]))),]
lmResults$dir <- as.numeric(as.character(lmResults$dir))
lmResults$p <- as.numeric(as.character(lmResults$p))
lmResults$adjP <- p.adjust(lmResults$p,method="fdr") # adjust all p-values for multiple hypothesis testing!
lmResults$mLog10P <- -log10(lmResults$adjP) * lmResults$dir # this is a transformation of the p-value so that smaller p-values get larger numbers (larger bars in the plot later) and so the negative values represent an inverse relationship and positive values a direct relationship between factor and effect
lmResults$color <- plotColLookup$datasetName[as.character(lmResults$dataset)] 
lmResults$variable <- recodeVarType(lmResults$variable)
lmResults$dataset <- factor(datasetTable[as.character(lmResults$dataset),"prettyLabel"],levels=(datasetTable[datasetNamesInOrder,"prettyLabel"]))

message("\t* make plots")

ggData <- lmResults[!is.na(lmResults$p),] #corResults[!is.na(corResults$p),] #

# the plot contains a few very extreme p-values, which need to be cropped so that the other p-values can be
# better visualized. We therefore trim mLog10P values above or below the 95% quantile (these values will be 
# annotated in the plot later on):
maxP <- max(-log10(0.05)*1.1,ceiling(quantile(abs(ggData$mLog10P),0.95,na.rm=T)+0.1))
ggData$mLog10P[ggData$mLog10P > maxP] <- maxP
ggData$mLog10P[ggData$mLog10P < -maxP] <- -maxP

# generate a panel of barplots:
d <- ggplot(ggData, aes(x=dataset,y=mLog10P,fill=color)) + geom_bar(width=0.5,colour="#333333",size=0.25,stat="identity")  + scale_fill_identity() + geom_hline(yintercept=0)  + facet_wrap( ~ variable, scales="free_y", ncol=1, drop=T) + defaultTheme(flipX=TRUE) + xlab(NULL)+ ylim(c(-maxP,maxP)) + geom_hline(yintercept=-log10(0.05),colour="gray",linetype=2) + geom_hline(yintercept=log10(0.05),colour="gray",linetype=2)

svgPlotGG(d, "dnaconc_influences", 11, 11)
# to find the results outside the upper lower transformed p-value bounds:
# (to be marked or corrected by editing the figure)
# lmResults[abs(lmResults$mLog10P)>maxP,]

# generate illustrative scatter plots for those assays that are most strongly influenced by each region characteristic:
ggScatterData <- melt(modelData[,c("datasetName",dependentVar,selColsNumeric)],id.vars=c("datasetName",dependentVar))
ggScatterData$datasetName <- factor(datasetTable[as.character(ggScatterData$datasetName),"prettyLabel"],levels=rev(datasetTable[datasetNamesInOrder,"prettyLabel"]))
for(curVar in selColsNumeric) {	
	tmp <- ggData[ggData$variable==recodeVarType(curVar),]
	for(curDataset in unique(tmp$dataset)) {
		tmp2 <- ggScatterData[ggScatterData$datasetName==curDataset & ggScatterData$variable==curVar,]
		maxY <- 100
		rangeX <- range(tmp2[,"value"], finite=T)*10
		rangeX <- c(0.9*(rangeX[1]),1.1*(rangeX[2]))/10

		svgPlot(paste0("dnaconc_scatter_influences_",curVar,"_",curDataset), 2.4, 2.6)
		par(cex=1)
		par(mar=c(4,4,2,1))
		plot(tmp2$value, abs(tmp2[,dependentVar]), pch=16, xlab=recodeVarType(curVar), ylab="Abs. dev. from corridor", las=1,col=rgb(0,0,0,0.3), ylim=c(0,maxY), xlim=rangeX, main=curDataset)
		void <- tryCatch(abline(lm(as.formula(paste0(dependentVar, " ~ value")), tmp2), col="hotpink"),error=function(e) {})
		par(xpd=NA)
		text(mean(rangeX),y=1.1 * maxY,sprintf("r = %.1f, adj.P = %.2f", cor(abs(tmp2[,dependentVar]), tmp2$value, use="pairwise"), tmp[tmp$dataset==curDataset,"adjP"]))
		dev.off()
	}
}





############# SPECIFIC NUMBERS AND STATISTICS FOR PAPER #############

message("Other statistics:")

# numbers quoted in section "Selection and provision of reference samples and regions"
message("\t* number of locus-specific measurements passing lab-QC: ",sum(!is.na(lsaData$methValue)),"\n")
message("\t* number of global measurements passing lab-QC: ",sum(!is.na(globalMethLevelsRaw)),"\n")

# numbers quoted in section "Performance of absolute assays"

absCorTable <- melt(corTable[datasetsByType$absolute,datasetsByType$absolute])
absCorTable <- absCorTable[as.numeric(absCorTable$Var1)<as.numeric(absCorTable$Var2),]
for(thresh in c(0.9,0.8)) {
	message("\t* number of absolute assay comparisons with correlation coefficient >",thresh,"= ",sum(absCorTable$value>thresh,na.rm=T) / nrow(absCorTable),"\n")
}
bestAbsCorAssays <- (absCorTable[absCorTable$value==max(absCorTable$value,na.rm=T),])
message("\t* test correlation: Assays=[",as.character(bestAbsCorAssays$Var1),as.character(bestAbsCorAssays$Var2),"], r =",round(bestAbsCorAssays$value,2),"\n")
for(queryAssays in list(c("Pyroseq_1","Pyroseq_2"),c("Pyroseq_1","EpiTyper_3"),c("Pyroseq_1","AmpliconBS_1"))) {
	message("\t* query correlation: Assays=[",queryAssays,"], r =",round(corTable[queryAssays[1],queryAssays[2]],2),"\n")
}

m100 <- mean(aggregate((min+max/2)~regionName,lsaData[lsaData$datasetName%in%datasetsByType$absolute & lsaData$sampleName=="Titration1_1_100",c("regionName","min","max")],mean)[,2])
m75 <- mean(aggregate((min+max/2)~regionName,lsaData[lsaData$datasetName%in%datasetsByType$absolute & lsaData$sampleName=="Titration1_2_75",c("regionName","min","max")],mean)[,2])
message("\t* consensus mean of all regions for 100% and 75% in Titration 1, m100 =",m100,", m75 =", m75,"\n")

devFromCorridor <- aggregate(abs(value)~datasetName,melt(lsaData[lsaData$datasetName%in%datasetsByType$absolute & lsaData$regionName%in%regionNames & lsaData$sampleName%in%setdiff(sampleNames,sampleNames[grep("Titration",sampleNames)]),][,c("datasetName","deviationCorridor")]),mean)
devFromCorridor <- devFromCorridor[order(devFromCorridor[,2]),]

biasFromCorridor <- aggregate(value~datasetName,melt(lsaData[lsaData$datasetName%in%datasetsByType$absolute & lsaData$regionName%in%regionNames & lsaData$sampleName%in%setdiff(sampleNames,sampleNames[grep("Titration",sampleNames)]),][,c("datasetName","deviationCorridor")]),mean)
biasFromCorridor <- biasFromCorridor[order((biasFromCorridor[,2])),]

# numbers quoted in section "Performance of relative assays"
relCorTable <- melt(corTable[datasetsByType$relative,datasetsByType$relative])
relCorTable <- relCorTable[as.numeric(relCorTable$Var1)<as.numeric(relCorTable$Var2),]
message("\t* correlation of relative assays, ranging from ",round(min(relCorTable$value),2)," to ",round(max(relCorTable$value),2),"\n")



### DATA FOR SUPPLEMENTARY TABLE S1 ###

# total number of valid region targets (at least one non-NA measurement) per assay:
designSuccessful <- aggregate(regionName~datasetName,lsaData[!is.na(lsaData$methValue),],function(x){length(unique(x))})
rownames(designSuccessful) <- designSuccessful[,1]
# total number of valid samples (at least one non-NA measurement) per assay:
samplesAttempted <- aggregate(sampleName~datasetName,lsaData[!is.na(lsaData$methValue),],function(x){length(unique(x))})
rownames(samplesAttempted) <- samplesAttempted[,1]

# total number of valid measurments per assay:
measurementsSuccessful <- aggregate(sampleName~datasetName,lsaData[!is.na(lsaData$methValue),],function(x){length(x)})
rownames(measurementsSuccessful) <- measurementsSuccessful[,1]
measurementsSuccessful["MS_HRM",2] <- measurementsSuccessful["MS_HRM",2]+18 # values "H" have been translated to NA's in the data table, however, these measurements have actually been successful!

lsaSummary <- data.frame("Dataset"=datasetNamesInOrder,"Amount of DNA"=datasetTable[datasetNamesInOrder,"amountOfDNA"],"SNP Detected"=datasetTable[datasetNamesInOrder,"snpDetected"],"Design Attempted"=datasetTable[datasetNamesInOrder,"designAttempted"],"Design Successful"=designSuccessful[datasetNamesInOrder,2],"Samples Attempted"=samplesAttempted[datasetNamesInOrder,2],"Measurements Attempted"=designSuccessful[datasetNamesInOrder,2]*samplesAttempted[datasetNamesInOrder,2],"Measurements Successful"=measurementsSuccessful[datasetNamesInOrder,2],check.names=FALSE)
lsaSummary <- cbind(
	lsaSummary,
	"Design Success Rate" = round(lsaSummary[,"Design Successful"]/lsaSummary[,"Design Attempted"]*100.0,1),
	"Assay Success Rate" = round(lsaSummary[,"Measurements Successful"]/lsaSummary[,"Measurements Attempted"]*100.0,1)
)

# general descriptive statistics about the measurements taken by all assays:
summaryStats <- list(
	"N"=function(x) { sum(!is.na(x)) }
	,"Mean"=function(x) { mean(x,na.rm=T) }
	,"10% Quantile"=function(x) { quantile(x,0.10,na.rm=T) }
	,"25% Quantile"=function(x) { quantile(x,0.25,na.rm=T) }
	,"Median / 50% Quantile"=function(x) { quantile(x,0.5,na.rm=T) }
	,"75% Quantile"=function(x) { quantile(x,0.75,na.rm=T) }
	,"90% Quantile"=function(x) { quantile(x,0.90,na.rm=T) }
	,"Standard Deviation"=function(x) { sd(x,na.rm=T) }
	,"Minimum"=function(x) { min(x,na.rm=T) }
	,"Maximum"=function(x) { max(x,na.rm=T) }
)
for(summaryStat in names(summaryStats)) {
	tmp <- aggregate(methValue~datasetName,lsaData, summaryStats[[summaryStat]] )
	rownames(tmp) <- tmp[,1]
	lsaSummary <- cbind(lsaSummary,tmp[datasetNamesInOrder,2])
	colnames(lsaSummary)[length(lsaSummary)] <- summaryStat
}

# descriptive statistics based on the deviation of measurements to the consensus:
summaryStatsConsensus <- list(
	"Deviation from Consensus - Abs. Mean"=function(x) { mean(abs(x),na.rm=T) }
	,"Deviation from Consensus - Abs. Median"=function(x) { median(abs(x),na.rm=T) }
	,"Deviation from Consensus - Bias"=function(x) { mean(x,na.rm=T) }
)
for(summaryStat in names(summaryStatsConsensus)) {
	tmp <- aggregate(deviationCorridor~datasetName,lsaData[lsaData$sampleName %in% setdiff(sampleNames,sampleNames[grep("Titration",sampleNames)]),], summaryStatsConsensus[[summaryStat]] )
	rownames(tmp) <- tmp[,1]
	lsaSummary <- cbind(lsaSummary,tmp[datasetNamesInOrder,2])
	colnames(lsaSummary)[length(lsaSummary)] <- summaryStat
}

# performance characteristics for titration series:
tmp <- performanceTableTitrations[
	regexpr("region_",performanceTableTitrations$regionName)>=0,
	c("datasetName","titrationType","Pearson.s.r","Adjusted.R.2","Residual.standard.error")
]
tmp <- melt(tmp,id.vars=c("datasetName","titrationType"))
tmp <- tmp[!is.na(tmp$value) & !is.nan(tmp$value),]
levels(tmp$variable)[levels(tmp$variable)=="Pearson.s.r"] <- "Pearson's r"
levels(tmp$variable)[levels(tmp$variable)=="Adjusted.R.2"] <- "Adjusted R^2"
levels(tmp$variable)[levels(tmp$variable)=="Residual.standard.error"] <- "Residual standard error"
for(titrationType in unique(tmp$titrationType)) {
	for(metricType in unique(tmp$variable)) {
		tmp2 <- aggregate(value~datasetName,tmp[tmp$titrationType==titrationType & tmp$variable==metricType,], function(x) { median(x,na.rm=T) } )
		rownames(tmp2) <- tmp2[,1]
		lsaSummary <- cbind(lsaSummary,tmp2[datasetNamesInOrder,2])
		colnames(lsaSummary)[length(lsaSummary)] <- paste(titrationType,"-",metricType,"(median)")
	} 
}

rownames(lsaSummary) <- lsaSummary[,1]

lsaSummary <- round(lsaSummary[,3:ncol(lsaSummary)],3)
lsaSummary[is.na(lsaSummary)] <- "n/a"

write.table(lsaSummary,file="results_analysis/misc/table-s1-stats.tsv",sep="\t",quote=FALSE,row.names=TRUE,col.names=NA)





