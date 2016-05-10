# SENSITIVITY ANALYSIS FOR CORRIDOR PARAMETERS

# to examine the effect of the choice of parameters for consensus corridors:
# 	- number of contributing technologies (default: CORRIDOR_MIN_NUM_ASSAY_TYPES)
# 	- width of flanking regions (default: +/- (CORRIDOR_EXTENSION/2))

# ensure reproducibility by fixing the random seed to an arbitrary offset:
load(latestPrepName)
set.seed(1234)
message("=== SENSITIVITY ANALYSIS ===")

# load library and define functions to be used:
library("sensitivity")
library("ggplot2")
library("pse")
library("reshape2")
library("parallel")

# set constants:
selDatasets <- datasetsByType$absolute
selRegions <- regionNames
selSamples <- setdiff(sampleNames,sampleNames[grep("Titration",sampleNames)])
originalLsaData <- lsaData[lsaData$datasetName%in%selDatasets & lsaData$regionName%in%selRegions & lsaData$sampleName%in%selSamples,]

# use the values from the actual benchmark study as a reference:
referenceValues <- aggregate(abs(deviationCorridor)~datasetName, originalLsaData, mean, na.rm=TRUE)
referenceValues <- structure(referenceValues[,2], names=referenceValues[,1])

# define parameter range to be examined:
parRanges <- list(
	minAssayTypeNumValues = 2:length(unique(datasetTable[datasetsByType$absolute, "assayGroup"])), # examine all combinations with at least two and at most all different technologies
	corridorExtValues = seq(0,50,by=5) # examine all possible corridor extension values from 0 to 50%, in 0.5% steps
)

### Define additional functions for sensitivity analysis: ###

# prepare measurement data using alternative corridor parameters:
prepData <- function(corridorMinAssayTypes=CORRIDOR_MIN_NUM_ASSAY_TYPES, corridorExtension=CORRIDOR_EXTENSION) {
	lsaData <- originalLsaData

	rNames <- unique(lsaData$regionName)
	sNames <- unique(lsaData$sampleName)

	consensusMeths <- list("lower"=function(x) { max(0,min(x,na.rm=T)-corridorExtension/2) },"upper"=function(x) { min(100,max(x,na.rm=T)+corridorExtension/2) },"median"=function(x) { median(x,na.rm=T) })

	for (curRegion in rNames) {
		curRegionData <- lsaData[lsaData$regionName==curRegion & lsaData$datasetName %in% datasetsByType$absolute,]
		curRegionData$assayGroup <- datasetTable[curRegionData$datasetName,"assayGroup"]
		for (curSample in sNames) {
			# subset the data for only those measurements concerning the current region/sample combination:
			curRegionSampleData <- curRegionData[curRegionData$sampleName==curSample,] 

			# figure out all possible intervals by looking up all values (+ the absolute lower and upper bound, -0.001 and 100),
			# sorting these values and then looking at all combinations:
			cutoffs <- sort(unique(c(-0.001,100,curRegionSampleData$methValue)))
			intervals <- data.frame(t(combn(cutoffs,2)))
			names(intervals) <- c("lower","upper")

			# sort the intervals by width and take only those that are wider than 0:
			intervals$width <- intervals$upper - intervals$lower
			intervals <- intervals[intervals$width>0,]
			intervals <- intervals[order(intervals$width),]

			# now iterate all intervals (starting with the smallest) until we find one containing at least 3 different assay types:
			for (i in 1:nrow(intervals)) {
				selData <- curRegionSampleData[curRegionSampleData$methValue>=intervals[i,"lower"] & curRegionSampleData$methValue<=intervals[i,"upper"],]
				if (length(unique(selData$assayGroup))>=corridorMinAssayTypes) break
			}

			# if we managed to find a valid interval, that's the consensus corridor:
			if (length(unique(selData$assayGroup))>=corridorMinAssayTypes) {	
				for (val in names(consensusMeths)) {
					lsaData[lsaData$regionName==curRegion & lsaData$sampleName==curSample,val] <- consensusMeths[[val]](selData$methValue)
				}
			}			
		}		
	}

	# calculate the deviation of each dataset's measurements from the borders of the consensus corridor:
	lsaData$deviationMedian <- lsaData$methValue - lsaData$median
	lsaData$deviationCorridor <- ifelse(lsaData$methValue>=lsaData$lower & lsaData$methValue<=lsaData$upper,0,ifelse(lsaData$methValue<lsaData$lower, lsaData$methValue-lsaData$lower, lsaData$methValue-lsaData$upper))

	return(lsaData)
}

# evaluate the results of a run using a variety of metrics:
evalRes <- function(x) {

	tmp <- aggregate(abs(deviationCorridor)~datasetName, x, mean, na.rm=TRUE)
	tmp <- sort(structure(tmp[,2],names=tmp[,1]))

	x <- as.numeric(tmp)	# the values according to the current analysis run
	y <- as.numeric(referenceValues[names(tmp)])	# the values according to the reference analysis run
	
	res <- list(
		"mean"=mean(x, na.rm=TRUE), # average value
		"pearson"=cor(x,y,method="pearson",use="pairwise"),
		"spearman"=cor(x,y,method="spearman",use="pairwise"),
		"best"=names(tmp)[1],
		"worst"=names(tmp)[length(tmp)]
	)

	for(deltaN in 1:floor(length(x)/2)) {
		top <- x[1:deltaN]
		bot <- x[(length(x)-deltaN+1):length(x)]
		res[[paste0("delta",deltaN)]] <- mean(top)-mean(bot)
	}

	append(res,tmp)
}

# run one analysis trial (examinining the performance of assays against the consensus corridor):
runAnalysisTrial <- function(corridorMinAssayTypes, corridorExtension) {
	message("corridorMinAssayTypes = ", corridorMinAssayTypes, ", corridorExtension = ", corridorExtension)

	curData <- prepData(corridorMinAssayTypes=corridorMinAssayTypes, corridorExtension=corridorExtension)
	curData <- curData[,c("datasetName","deviationCorridor")]

	evalRes(curData)
}





### Run the actual analysis and generate plots ###

## 1. Parameter space exploration with defined values:

pseResults <- list()
# loop through corridor parameter values to be explored and evaluate the performance:
for(minAssayTypeNum in parRanges$minAssayTypeNumValues) {
	pseResults[[as.character(minAssayTypeNum)]] <- list()
	for(corridorExt in parRanges$corridorExtValues) {
		message("corridorMinAssayTypes = ", minAssayTypeNum, ", corridorExtension = ", corridorExt)

	
		curData <- prepData(corridorMinAssayTypes=minAssayTypeNum, corridorExtension=corridorExt)
		e <- evalRes(curData[,c("datasetName","deviationCorridor")])

		pseResults[[as.character(minAssayTypeNum)]][[as.character(corridorExt)]] <- e
	}
} 

# prepare plots:
ggData <- melt(pseResults)
colnames(ggData) <- c("value","metric","corridorExtValues","minAssayTypeNumValues")
ggData$corridorExtValues <- factor(as.numeric(ggData$corridorExtValues)/2, levels=parRanges$corridorExtValues/2)
ggData$minAssayTypeNumValues <- factor(ggData$minAssayTypeNumValues, levels=parRanges$minAssayTypeNumValues)

# ... line plots illustrating assay rankings after evaluation with different corridors:
ggData2 <- ggData[ggData$metric%in%datasetNames,]
ggData2$value <- as.numeric(ggData2$value)
ggData2$corridorExtValues <- as.numeric(as.character(ggData2$corridorExtValues))
ggData2$minAssayTypeNumValues <- paste("Min. assay types: ", ggData2$minAssayTypeNumValues)
ggData2$datasetName <- factor(ggData2$metric, levels=names(sort(referenceValues)), labels=datasetTable[names(sort(referenceValues)),"prettyLabel"])
refDispData <- data.frame(corridorExtValues=CORRIDOR_EXTENSION,minAssayTypeNumValues=paste("Min. assay types: ", CORRIDOR_MIN_NUM_ASSAY_TYPES),datasetName=datasetTable[names(referenceValues),"prettyLabel"],value=referenceValues)
d <- ggplot(ggData2, aes(datasetName, value, group=corridorExtValues, color=corridorExtValues)) + facet_grid(~minAssayTypeNumValues) + geom_line() + defaultPlotTheme(flipX=TRUE, fontSize=8) + annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + xlab("Dataset (ordered by average absolute deviation with default corridor parameters)") + ylab("Average abs. deviation from corridor") + scale_color_gradient2(low="#ffeda0",mid="#31a354",high="#002200",midpoint=12.5) + geom_line(data=refDispData, color="red")
svgPlotGG(d, "sensitivity_2b_lines", 22, 8, units="cm")

# ... best/worst assay per trial run:
ggData2 <- ggData[ggData$metric%in%c("best","worst"),]
ggData2$assayGroup <- datasetTable[as.character(ggData2$value),"assayGroup"]
ggData2$prettyLabel <- datasetTable[as.character(ggData2$value),"prettyLabel"]
ggData2$metric <- factor(ggData2$metric, levels=c("best","worst"), labels=c("Best (lowest average abs. deviation)", "Worst (highest average abs. deviation)"))
d <- ggplot(ggData2, aes(corridorExtValues, minAssayTypeNumValues)) + geom_tile(aes(fill = assayGroup), colour = "white") + scale_fill_manual(values=plotColLookup$assayGroup, guide=FALSE) + geom_text(aes(label=prettyLabel),family="Arial",colour="black",size=2.6) + ylab("Min. number of assay types in corridor") + xlab("Corridor extension (+/- X%)") + coord_flip()  + facet_grid(metric~., scales="free", margins=F) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),	panel.margin = unit(0, "lines"), axis.title= element_text(size=8,face="bold"),	axis.text.y = element_text(hjust=1,vjust=0.5,size=8), axis.text.x = element_text(hjust=1,vjust=0.5,size=8), axis.ticks = element_blank(), strip.background = element_blank(), strip.text = element_text(size=8,face="bold"))
svgPlotGG(d, "sensitivity_2b_bestworst", 11, 14, units="cm")





## 2. Latin Hypercube Sampling:

# perform sampling runs:
myLHS <- LHS(model=function(X) { print(X); return(mcmapply(function(a,b) { runAnalysisTrial(a,b) }, round(X[,1]), X[,2], mc.cores=16)) }, factors=names(parRanges), N=1000, q=c("qunif","qunif"), q.arg=lapply(parRanges,function(x) list("min"=min(x),"max"=max(x))), nboot=50)

# do one trial run to get names:
defRun <- runAnalysisTrial(CORRIDOR_MIN_NUM_ASSAY_TYPES,CORRIDOR_EXTENSION)
resNames <- names(defRun)

# prepare plots:
ggData <- data.frame(get.data(myLHS),as.data.frame(get.results(myLHS,get.mean=FALSE)))
colnames(ggData) <- c("Min. number of assay types in corridor", "Corridor extension (+/- X%)", resNames)
ggData[,1] <- round(ggData[,1])
ggData[,2] <- ggData[,2]/2
ggData[,-(1:2)] <- apply(ggData[,-(1:2)],2,unlist)
ggData <- melt(ggData, id.vars=c("Min. number of assay types in corridor","Corridor extension (+/- X%)"))

ggData2 <- ggData[ggData$variable%in%c("mean","pearson","spearman"),]
ggData2$value <- as.numeric(ggData2$value)
ggData2$variable <- factor(ggData2$variable, levels=c("mean","pearson","spearman"), labels=c("Average abs. deviation","Pearson correlation vs. default","Spearman correlation vs. default"))
ylims <- c(10,1,1)

# plot evaluation by all available metrics:
for(i in 1:length(levels(ggData2$variable))) {
	# performance by minimum number of contributing assays:
	d <- ggplot(ggData2[ggData2$variable==levels(ggData2$variable)[i],], aes(`Min. number of assay types in corridor`, value, color=`Corridor extension (+/- X%)`)) + geom_point(size=0.8) + ylab(levels(ggData2$variable)[i]) + defaultPlotTheme() + scale_color_gradient2(low="#ffeda0",mid="#31a354",high="#002200",midpoint=mean(ggData$`Corridor extension (+/- X%)`),guide=guide_colorbar(title="Colors:",direction="horizontal")) + ylim(0,ylims[i]) + theme(panel.grid.major=element_line(color="lightgray",linetype=2),panel.grid.major.x=element_blank()) + geom_point(x=CORRIDOR_MIN_NUM_ASSAY_TYPES,y=as.numeric(defRun[i]),color="red",size=4,pch=4)
	svgPlotGG(d, paste0("sensitivity_lhs_byminassay_",i), 10, 7, units="cm")

	# performance by corridor extension:
	d <- ggplot(ggData2[ggData2$variable==levels(ggData2$variable)[i],], aes(`Corridor extension (+/- X%)`, value, color=as.factor(`Min. number of assay types in corridor`))) + geom_point(size=0.8) + ylim(0,ylims[i]) + ylab(levels(ggData2$variable)[i]) + defaultPlotTheme() + scale_color_manual(values=c("#018571","#AAD6CF","#EFB3D8","#D01C8B"),guide=guide_legend(title="Colors:",direction="horizontal")) + theme(panel.grid.major=element_line(color="lightgray",linetype=2),panel.grid.major.x=element_blank()) + geom_point(x=CORRIDOR_EXTENSION/2,y=as.numeric(defRun[i]),color="red",size=4,pch=4)
	svgPlotGG(d, paste0("sensitivity_lhs_bycorridorext_",i), 12, 7, units="cm")
}


















