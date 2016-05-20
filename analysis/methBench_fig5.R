# FIGURE 5: PERFORMANCE OF LOCAL DNA METHYLATION ASSAYS FOR TRAINING CLASSIFIERS DISCRIMINATING NORMAL FROM TUMOR TISSUE


# ensure reproducibility by fixing the random seed to an arbitrary offset:
load(latestPrepName)
set.seed(1234)
message("=== FIGURE 5 ===")




# define comparisons to be examined:
cmp <- list("FG"=c(sampleNames[grepl("Tumor",sampleTable$sampleName) & sampleTable$sampleType=="TumorNormal"]),"BG"=c(sampleNames[grepl("Normal",sampleTable$sampleName) & sampleTable$sampleType=="TumorNormal"]))
combos <- numeric()
for(testIndex1 in 1:(length(cmp$FG)-1)) {
	for(testIndex2 in (testIndex1+1):length(cmp$FG)) {
		combos <- rbind(combos,c(cmp$FG[c(testIndex1,testIndex2)],cmp$BG[c(testIndex1,testIndex2)]))
	}
}


# choose data to work with:
selDatasets <- datasetNames
selRegions <- coreRegionNames # focus on core regions, that have been measured by almost all assays!






### ROC CURVES FOR NOISE/ERROR MODELS (Figure 5A) ###

message("Train tumor/normal classifiers with noise...")
rocPlotDataClassifiers <- list()

# train support vector classifiers under various noise/error models and levels based on the data 
# of each individual local DNA methylation assay. For training/testing, use a 2-fold cross-validation
# procedure, in which two pairs of sample (2x matched tumor/normal tissue samples) are held out for
# model evaluation and the rest is used for training: 
# (plot ROC curves as we go along!)
performanceTableClassifiers <- numeric()
for(curDataset in selDatasets) {	
	# print progress info:	
	message("\t* ", Sys.time(), " ", curDataset," (",which(curDataset==selDatasets)," of ",length(selDatasets),")")

	# find fore- and background data for the current dataset:
	fgI <- which(mlData$dataset==curDataset & mlData$sample%in%cmp$FG)
	bgI <- which(mlData$dataset==curDataset & mlData$sample%in%cmp$BG)

	# put both sets of data together into one table:
	mlDataX <- rbind(mlData[fgI,selRegions],mlData[bgI,selRegions])
	rownames(mlDataX) <- mlData$sample[c(fgI,bgI)]
	mlDataX <- mlDataX[c(cmp$FG,cmp$BG),]
	rownames(mlDataX) <- c(cmp$FG,cmp$BG)
	mlDataX <- mlDataX[,colSums(!is.na(mlDataX))>1] # only keep regions with at least 2 measurements

	# extend the table by prepending a column with the target labels, i.e. tumor or normal:
	mlDataX <- cbind("Target"=c(rep(1,length(cmp$FG)),rep(0,length(cmp$BG))),mlDataX)
	mlDataX <- as.data.frame(mlDataX)

	# plot ROC curves with model performance results:
	svgPlot(paste0("roc_",curDataset), 6, 6)
	first <- TRUE
	rocPlotDataTmp <- list()

	# then train classifiers with the different noise/error models...
	for(noiseType in names(NOISE_TYPES)) {	

		noiseFun <- NOISE_TYPES[[noiseType]]
		rocPlotDataTmpN <- list()	

		# ... and at different noise/error levels:
		for(noiseLevelI in 1:length(NOISE_LEVELS)) {
			noiseLevel <- NOISE_LEVELS[noiseLevelI]

			# avoid doing additional work, so only train classifiers at noise/error level=0 once:
			if(noiseLevel==0 & noiseType!=names(NOISE_TYPES)[1]) next;

			# print progress info:
			lbl <- paste(curDataset,', n=',noiseLevel,sep="")
			message("\t\t",noiseType,"\t",lbl)

			# run classifier trials:
			trialResult <- runClassifierTrials(mlDataX,rndDataModFun=function(x){noiseFun(x,n=noiseLevel)},nTrials=CLASSIFIER_TRIALS)

			# record performance metrics (AUC = area under curve):
			performanceTableClassifiers <- rbind(performanceTableClassifiers,c(
				"dataset" = curDataset,
				"noise" = paste(noiseType,noiseLevel,collapse=":"),
				"noiseLevel" = noiseLevel,
				"noiseType" = noiseType, 
				"AUC"=trialResult$auc
			))

			# at a line for the current noise level to the ROC plot:
			plot(trialResult$perf, col=NOISE_COLORS[[noiseType]][noiseLevelI],add=!first, lty=which(NOISE_LEVELS==noiseLevel), main=datasetTable[curDataset,"prettyLabel"],lwd=2.75, avg="threshold")
			first <- F

			rocPlotDataTmpN[[noiseLevelI]] <- trialResult$perf
		}
		names(rocPlotDataTmpN) <- as.character(NOISE_LEVELS)

		rocPlotDataTmp[[noiseType]] <- rocPlotDataTmpN
	}

	# add legends for colors and line types to the ROC plot:
	legend("bottomright",as.character(NOISE_LEVELS),lty=1:length(NOISE_LEVELS), bty = "n", cex = 1, title="Noise Level",lwd=2.75)
	cols <- c()
	for(noiseType in names(NOISE_TYPES)) {
		cols <- c(cols,NOISE_COLORS[[noiseType]][2])
	}
	legend("bottom",as.character(names(NOISE_TYPES)),col=cols, cex = 1, title="Noise Type", bty = "n", lwd=2.75)

	rocPlotDataClassifiers[[curDataset]] <- rocPlotDataTmp
	dev.off()
}

# find best predictor (mean of AUCS --> robust to noise):
aucData <- as.data.frame(performanceTableClassifiers)[performanceTableClassifiers[,"noiseType"] %in% c("Random error","Uniform noise"),]
aucData[,"AUC"] <- as.numeric(aucData[,"AUC"])
aucsNoise <- aggregate(AUC~dataset,data=aucData,mean,na.rm=T)
aucsNoise <- aucsNoise[order(aucsNoise$AUC),]


### SUPPLEMENTARY PANEL OF ROC CURVES (Supplementary Figure 14) ###

svgPlot("roc_small_all", 2*3,2.2*7, pointsize=18)
par(mfrow=c(7,3))
par(mar=c(3,2,1,1))
par(mgp=c(3,0.75,0))
for(curDataset in selDatasets) {
	plotSelectedNoiseRocs(curDataset, withLegend=which(curDataset==selDatasets)==length(selDatasets))
}
dev.off()




### AUC SUMMARY SUBJECT TO NOISE/ERROR MODELS (Figure 5B) ###

message("AUC summary plots for tumor/normal classifiers with noise (Figure 5.b)... ")

# for classifier training with noise/error level =0, we have only trained a model for "Random error", but not for "Uniform noise",
# since both use the original, unperturbed data in this case. For plotting purposes, we need to replicate the 0-error results for 0-noise,
# so we extend the original data table with a copy of the respective results:
tmp <- performanceTableClassifiers[performanceTableClassifiers[,"noiseLevel"]==0&performanceTableClassifiers[,"noiseType"]=="Random error",]
tmp[,"noiseType"] <- "Uniform noise"
tmp <- rbind(performanceTableClassifiers,tmp)

# then plot the evaluation results:
svgPlot("auc-change", 10, 6)
plotAUCChange(tmp,by="noiseLevel",categoryTypesName="noiseType",selectedDatasets=selDatasets,exclude=c(1),asc=T)
dev.off()




### ROC CURVES FOR FEATURE SELECTION MODELS (Figure 5C) ###

message("Train tumor/normal classifiers with selected features...")
rocPlotDataClassifiersSelection <- list()

# like before, we train SVM classifiers in a cross-validation procedure, but instead of examining different noise/error models,
# we now select a subset of features (= regions) to train and evaluate the models with. The features are either chosen randomly
# or by examining the likely information content using an F-score.

noiseFun <- NOISE_TYPES[[STANDARD_NOISE_TYPE]] # keep adding a small amount of noise also for this experiment
selectionColors <- brewer.pal(3,"PiYG")[c(1,3)]
performanceTableClassifiersSelection <- numeric()
for(curDataset in selDatasets) { 	
	# print progress info:	
	message("\t* ", Sys.time(), " ", curDataset," (",which(curDataset==selDatasets)," of ",length(selDatasets),")")

	# find fore- and background data for the current dataset:
	fgI <- which(mlData$dataset==curDataset & mlData$sample%in%cmp$FG)
	bgI <- which(mlData$dataset==curDataset & mlData$sample%in%cmp$BG)

	# put both sets of data together into one table:
	mlDataX <- rbind(mlData[fgI,selRegions],mlData[bgI,selRegions])
	rownames(mlDataX) <- mlData$sample[c(fgI,bgI)]
	mlDataX <- mlDataX[c(cmp$FG,cmp$BG),]
	rownames(mlDataX) <- c(cmp$FG,cmp$BG)
	mlDataX <- mlDataX[,colSums(!is.na(mlDataX))>1] # only keep regions with at least 2 measurements

	# extend the table by prepending a column with the target labels, i.e. tumor or normal:
	mlDataX <- cbind("Target"=c(rep(1,length(cmp$FG)),rep(0,length(cmp$BG))),mlDataX)
	mlDataX <- as.data.frame(mlDataX)

	# plot ROC curves with model performance results:
	svgPlot(paste0("roc_selection_",curDataset), 6, 6)
	first <- TRUE
	rocPlotDataTmp <- list()

	# train classifiers with the different feature selection methods...
	for(regionSelectionMethodName in names(REGION_SELECTION_METHODS)) {
		rocPlotDataTmpN <- list()	

		# ... and using different numbers of features (regions):
		for(howMany in REGION_NUMBERS) {
			regionSelectionMethod <- function(x) { noiseFun(REGION_SELECTION_METHODS[[regionSelectionMethodName]](x,howMany=howMany),n=STANDARD_NOISE) }
	
			# print progress info:				
			lbl <- paste(curDataset,', n=',regionSelectionMethodName,"-",howMany,sep="")
			message("\t\t",regionSelectionMethodName,"-",howMany)

			# run classifier trials:
			trialResult <- runClassifierTrials(mlDataX,rndDataModFun=regionSelectionMethod,nTrials=CLASSIFIER_TRIALS)

			# record performance metrics (AUC = area under curve):
			performanceTableClassifiersSelection <- rbind(performanceTableClassifiersSelection,c(
				"dataset" = curDataset,
				"selectionFun" = regionSelectionMethodName,
				"regNumber" = howMany,
				"AUC"=trialResult$auc
			))

			# at a line for the current feature selection method to the ROC plot:
			plot(trialResult$perf, col=selectionColors[which(regionSelectionMethodName==names(REGION_SELECTION_METHODS))],add=!first, lty=which(howMany==REGION_NUMBERS), main=datasetTable[curDataset,"prettyLabel"],lwd=2.75, avg="threshold")	
			first <- FALSE

			rocPlotDataTmpN[[which(howMany==REGION_NUMBERS)]] <- trialResult$perf
		}
		names(rocPlotDataTmpN) <- as.character(REGION_NUMBERS)

		rocPlotDataTmp[[regionSelectionMethodName]] <- rocPlotDataTmpN
	}

	# add legends for colors and line types to the ROC plot:
	legend("bottomright",as.character(REGION_NUMBERS),lty=1:length(REGION_NUMBERS), bty = "n", cex = 1, title="Number of regions",lwd=2.75)
	legend("bottom",as.character(names(REGION_SELECTION_METHODS)),col=selectionColors, cex = 1, title="Selection method", bty = "n", lwd=2.75)
	dev.off()

	rocPlotDataClassifiersSelection[[curDataset]] <- rocPlotDataTmp

}


### SUPPLEMENTARY PANEL OF ROC CURVES (Supplementary Figure 14) ###

svgPlot("roc_small_all_selection", 2*3,2.2*7, pointsize=18)
par(mfrow=c(7,3))
par(mar=c(3,2,1,1))
par(mgp=c(3,0.75,0))
for(curDataset in selDatasets) {
	plotSelectedSelectionRocs(curDataset, withLegend=which(curDataset==selDatasets)==length(selDatasets))
}
dev.off()


### AUC SUMMARY SUBJECT TO FEATURE SELECTION METHODS (Figure 5D) ###

message("AUC summary plots for tumor/normal classifiers with selected features (Figure 5.d)... ")

svgPlot("auc-change_selection", 10, 6)
plotAUCChange(performanceTableClassifiersSelection,by="regNumber",categoryTypesName="selectionFun",categoryTypes=c("Top","Random"),selectedDatasets=selDatasets,exclude=c())
dev.off()

# find best predictor (mean of AUCS --> robust to noise):
tmp <- as.data.frame(performanceTableClassifiersSelection)
tmp$AUC <- as.numeric(as.vector(tmp$AUC))
aucsSelection <- aggregate(AUC~dataset,data=tmp,mean,na.rm=T)
aucsSelection <- aucsSelection[order(aucsSelection$AUC),]
aucsSelection



save.image(paste("methBench","afterFig5",format(Sys.time(), "%Y%m%d_%H%M"),"session.Rdata",sep="_"))
