# FIGURE 4: GLOBAL DNA METHYLATION ANALYSIS AND PREDICTION


# ensure reproducibility by fixing the random seed to an arbitrary offset:
load(latestPrepName)
set.seed(1234)
message("=== FIGURE 4 ===")

# choose data to work with:
selSamples <- sampleNames
selRegions <- coreRegionNames
selDatasets <- datasetsByType$absoluteAndRelative




### GLOBAL DNA METHYLATION MEASUREMENTS (Figure 4A) ###

message("Barplots for global DNA methylation assays (Figure 4.a)...")
ggData <- melt(globalMethLevelsRaw)
ggData$color <- plotColLookup$sampleName[as.character(ggData$Var1)]
d <- ggplot(ggData, aes(x=Var1,y=value,fill=color)) + geom_bar(width=0.5,colour="#333333",size=0.25,stat="identity") + geom_hline(yintercept=0) + scale_fill_identity()  + facet_grid(Var2 ~ ., scales="fixed", margins=F) + ylab(DNA_METH_LABEL) + xlab(NULL) + defaultPlotTheme(flipX=TRUE) + theme(axis.line.x=element_blank())
svgPlotGG(d,"values_global", 9, 12)




### ESTIMATE GLOBAL DNA METHYLATION FROM LOCAL VALUES (for Figures 4B-4E) ###

message("Build models to infer global methylation from local measurements...")

# define all simple calculation methods to be evaluated:
# (each method takes a set of local DNA methylation measurements as input and 
#  calculcates an estimate of the global DNA methylation level from those)
estimatorFuns <- list("mean"=mean,"median"=median)

# define the target values that should be predicted (for training regression models):
# 1. use an artifically crated "global target" reference as the average of the two best-performing assays
#    with some additional tweaks (see makeGlobalReference-method)
targetVals <- list("Global target" = makeGlobalReference(globalMethLevelsRaw,c("HPLC/MS","Pyroseq NBL2")))
# 2. also use each of the global DNA methylation assays individually as a reference
for(refName in colnames(globalMethLevelsRaw)) {
	targetVals[[refName]] <- globalMethLevelsRaw[mlData$sample,refName]
}

# now train regression models and evaluate predictions:
predictorResults <- numeric(0)
allPredictions <- matrix(ncol=3+length(selSamples),nrow=0,dimnames=list(NULL,c("predictor","dataset","targetType",as.character(selSamples))))
for(curDataset in selDatasets) { 

	message("\t* ",curDataset)

	# get the data for the current dataset:
	mlDataX <- mlData[mlData$dataset%in%ifelse(curDataset=="all",selDatasets,c(curDataset)),c("sample",selRegions)]
	mlDataX <- mlDataX[,colSums(!is.na(mlDataX))>1] # discard columns without valid measurements
	mlDataX[,2:ncol(mlDataX)] <- impute(mlDataX[,2:ncol(mlDataX)],what="median") # impute single missing measurements from the other valid measurements
	mlDataX <- as.data.frame(mlDataX)

	message("\t\t observations (sample/region pairs) = ",nrow(mlDataX),", features (regions) = ",ncol(mlDataX)-1)

	# then train models for each of the tested target values separately:
	for(targetValFunName in names(targetVals)) {
		message("\t\t ",targetValFunName)

		# build a combined data frame of target values and training data:
		mlDataComplete <- as.data.frame(cbind("Target"=targetVals[[targetValFunName]][as.character(mlDataX[,1])], data.matrix(mlDataX[,2:ncol(mlDataX)])))
		
		# ... and train and evaluate various types of simple estimator functions (mean, median) and
		# regression models (GLM, SVM, RF) using leave-one-out-cross-validation (LOOCV):

		# SIMPLE ESTIMATOR CALCULATIONS
		for(estimatorFunName in names(estimatorFuns)) {
			estimatorFun <- estimatorFuns[[estimatorFunName]]

			# calculcate predicted global DNA methylation values:
			pred <- apply(mlDataComplete[,2:ncol(mlDataComplete)],1,estimatorFun)
			vs <- structure(rep(NA,length(selSamples)),names=selSamples)
			vs[rownames(mlDataComplete)] <- pred

			# record predictions:
			allPredictions <- rbind(allPredictions,c("predictor"=estimatorFunName,"dataset"=curDataset,"targetType"=targetValFunName,vs))

			# evaluate predicted values vs. target values:
			evalResults <- evaluateRegressionResults(mlDataComplete$Target, pred)

			# record evaluation results:
			for(metricType in names(evalResults)) {
				predictorResults <- rbind(predictorResults, cbind("predictor"=estimatorFunName,"dataset"=curDataset,"targetType"=targetValFunName,"metricType"=metricType,evalResults[metricType]))
			}
		}

		# GENERALIZED LINEAR MODEL
		for(glm.family in c("gaussian")) {
			# train model and calculcate predicted global DNA methylation values:
			preds <- loocv(learner("runGLM",pars=list(fam=glm.family)), dataset(Target ~ ., mlDataComplete), loocvSettings())
			preds <- attr(preds, "foldResults")
			vs <- structure(rep(NA,length(selSamples)),names=selSamples)
			vs[rownames(mlDataComplete)] <- preds[,2]

			# record predictions:
			predictorName <- paste("GLM (",glm.family,")",sep="",collapse="")
			allPredictions <- rbind(allPredictions,c("predictor"=predictorName,"dataset"=curDataset,"targetType"=targetValFunName,vs))

			# evaluate predicted values vs. target values:
			evalResults <- evaluateRegressionResults(preds[,1],preds[,2])

			# record evaluation results:
			for(metricType in names(evalResults)) {
				predictorResults <- rbind(predictorResults, cbind("predictor"=predictorName,"dataset"=curDataset,"targetType"=targetValFunName,"metricType"=metricType,evalResults[metricType]))
			}
		}
	
		# SUPPORT VECTOR REGRESSION
		for(svm.kernel in c("linear", "polynomial")) {
			# train model and calculcate predicted global DNA methylation values:
			preds <- loocv(learner("runSVM",pars=list(k=svm.kernel)), dataset(Target ~ ., mlDataComplete), loocvSettings())
			preds <- attr(preds, "foldResults")
			vs <- structure(rep(NA,length(selSamples)),names=selSamples)
			vs[rownames(mlDataComplete)] <- preds[,2]

			# record predictions:
			predictorName <- paste("SVM (",svm.kernel,")",sep="",collapse="")
			allPredictions <- rbind(allPredictions,c("predictor"=predictorName,"dataset"=curDataset,"targetType"=targetValFunName,vs))

			# evaluate predicted values vs. target values:
			evalResults <- evaluateRegressionResults(preds[,1],preds[,2])

			# record evaluation results:
			for(metricType in names(evalResults)) {
				predictorResults <- rbind(predictorResults, cbind("predictor"=predictorName,"dataset"=curDataset,"targetType"=targetValFunName,"metricType"=metricType,evalResults[metricType]))
			}
		}

		# RANDOM FOREST REGRESSION

		# train model and calculcate predicted global DNA methylation values:
		preds <- loocv(learner("runRF",pars=list()), dataset(Target ~ ., mlDataComplete), loocvSettings())
		preds <- attr(preds, "foldResults")
		vs <- structure(rep(NA,length(selSamples)),names=selSamples)
		vs[rownames(mlDataComplete)] <- preds[,2]

		# record predictions:
		predictorName <- "Random Forest"
		allPredictions <- rbind(allPredictions,c("predictor"=predictorName,"dataset"=curDataset,"targetType"=targetValFunName,vs))

		# evaluate predicted values vs. target values:
		evalResults <- evaluateRegressionResults(preds[,1],preds[,2])

		# record evaluation results:
		for(metricType in names(evalResults)) {
			predictorResults <- rbind(predictorResults, cbind("predictor"=predictorName,"dataset"=curDataset,"targetType"=targetValFunName,"metricType"=metricType,evalResults[metricType]))
		}
	}
}
colnames(predictorResults) <- c("predictor","dataset","targetType","metricType","value")




### PERFORMANCE SUMMARY OF CLASSIFIERS (Figure 4C) ###

message("Classifier performance summary (Figure 4.c)...")

# prepare data for plotting:
dfPredictorResults <- as.data.frame(predictorResults,row.names=NA)
dfPredictorResults$value <- as.numeric(dfPredictorResults$value)
dfPredictorResults <- cbind(dfPredictorResults, color=as.character(plotColLookup$datasetName[dfPredictorResults$dataset]))
dfPredictorResults$dataset[dfPredictorResults$dataset %in% datasetsByType$absoluteAndRelative] <- datasetTable[dfPredictorResults$dataset[dfPredictorResults$dataset %in% datasetsByType$absoluteAndRelative],"prettyLabel"]
dfPredictorResults$dataset <- factor(dfPredictorResults$dataset,datasetTable[datasetNamesInOrder,"prettyLabel"])
selTTypes <- c("Global target",colnames(globalMethLevelsRaw))
dfPredictorResultsFiltered <- dfPredictorResults[dfPredictorResults$dataset!="all" & dfPredictorResults$predictor!="LM" & dfPredictorResults$targetType %in% selTTypes & dfPredictorResults$metricType %in% c("rmse"),]
dfPredictorResultsFiltered$metricType <- recode(dfPredictorResultsFiltered$metricType,"'r2'='R^2';'mae'='Mean Absolute Error';'rmse'='Root Mean Square Error'")
dfPredictorResultsFiltered$targetType <- factor(dfPredictorResultsFiltered$targetType,levels=selTTypes)
dfPredictorResultsFiltered$predictor <- factor(dfPredictorResultsFiltered$predictor,levels=unique(dfPredictorResultsFiltered$predictor))

# plot a panel of all prediction performance metrics (regression performance vs. global target and vs. individual global assays -
# this plot is not used in the paper, only for general interest):
d <- ggplot(dfPredictorResultsFiltered, aes(x=predictor,y=value)) + geom_boxplot() + geom_hline(y_intercept=0)  + facet_grid(targetType ~ metricType, scales="free")  + defaultPlotTheme(flipX=TRUE) + theme(axis.line.x=element_blank())
svgPlotGG(d, "globalMeth_prediction_small", 5, 22)
bestPredictors <- levels(dfPredictorResultsFiltered$predictor)[aggregate(value~targetType,aggregate(value~predictor+targetType,dfPredictorResultsFiltered,median),function(x) {y<-min(x); which(x==y)})$value]

# plot the performance metrics vs. the global target values (this is Fig. 4C):
d <- ggplot(dfPredictorResultsFiltered[dfPredictorResultsFiltered$targetType=="Global target",], aes(x=predictor,y=value)) + geom_boxplot() + ylab("Root mean square error") + xlab(NULL)  + defaultPlotTheme(flipX=TRUE)
svgPlotGG(d, "globalMeth_prediction_small_globalTarget", 5, 7)




### PERFORMANCE OF RANDOMFOREST REGRESSION WITH DIFFERENT DATASETS (Fig. 4D) ###

message("Performance of best model with different assays (Figure 4.d)...")
d <- ggplot(dfPredictorResultsFiltered[dfPredictorResultsFiltered$targetType=="Global target"&dfPredictorResultsFiltered$predictor=="Random Forest",], aes(x=dataset,y=value,fill=color)) + geom_bar(stat="identity",colour="#333333",width=0.65) + scale_fill_identity() + ylab("Root mean square error") + xlab(NULL) + defaultPlotTheme(flipX=TRUE)  + geom_hline(y_intercept=0) + theme(axis.line.x=element_blank())
svgPlotGG(d, "globalMeth_prediction_small_globalTarget_byAssay", 12, 8)




### CORRELATION OF PREDICTED AND REAL GLOBAL DNA METHYLATION MEASUREMENTS (Fig. 4B) ###

message("Correlation matrix: predictions and real measurements (Figure 4.b)...")
dfPredictions <- data.frame(allPredictions)
dfPredictions <- dfPredictions[dfPredictions$targetType=="Global target",]
predictors <- unique(dfPredictions$predictor)
for(curPredictor in predictors) {
	for(curTarget in unique(dfPredictions$targetType)) {
		curPreds <- dfPredictions[dfPredictions$predictor==curPredictor & dfPredictions$targetType==curTarget,]

		# extend the global DNA methylation data table with columns for the predictions derived from the current regression model 
		# using the current prediction method:
		globalMethExtended <- cbind(globalMethLevelsRaw,matrix(NA,nrow=nrow(globalMethLevelsRaw),ncol=length(selDatasets), dimnames=list(rownames(globalMethLevelsRaw),selDatasets)))
		for(curDataset in selDatasets) {
			nr <- nrow(curPreds[curPreds$dataset==curDataset,selSamples])
			if(nr>1) cat(curPredictor,curTarget,curDataset,nr,"\n")
			tmp <- colMeans(data.matrix(curPreds[curPreds$dataset==curDataset,selSamples]))
			globalMethExtended[selSamples,curDataset] <- tmp
		}		
		colnames(globalMethExtended)[colnames(globalMethExtended)%in%rownames(datasetTable)] <- datasetTable[colnames(globalMethExtended)[colnames(globalMethExtended)%in%rownames(datasetTable)],"prettyLabel"]
		
		# extend the data table with a column for the "global target" values":		
		globalMethExtended <- cbind("Global target"=
makeGlobalReference(globalMethLevelsRaw,c("HPLC/MS","Pyroseq NBL2")),globalMethExtended)

		
		# scale the values to the interval [0,100]
		globalMethExtended <- apply(globalMethExtended,2,scaleByMax)

		# calculcate the correlation table:
		# (N.B. prior to calculating the correlation coefficients, we impute missing values
		#  by finding the most closely matching other dataset with values and copying those values
		#  across. Datasets included in the calculation for the current correlation coefficient are
		#  excluded as possible imputation sources. We found that this way of dealing with missing
		#  values removed (positive and negative) biases affecting assays that had not attempted 
		#  all regions/samples)
		globalMethCor <- nnImputedCor(globalMethExtended,method="pearson")

		# plot the correlation table:
		svgPlot(paste0("cor_heatmap_predictor_estimates_nnImpute_",curPredictor,"_",gsub(":","",gsub("/","-",curTarget))), nrow(globalMethCor)/2, nrow(globalMethCor)/2)
		pheatmap(globalMethCor,col=colorRampPalette(rev(brewer.pal(5,"Oranges")),bias=1,space="Lab")(22),border_color="white",cluster_cols=TRUE,cluster_rows=TRUE,display_numbers=TRUE,cellwidth=22,cellheight=22,fontsize=13,breaks=seq(-0.1,1,0.05))
		dev.off()
	}
}




### PREDICTED GLOBAL DNA METHYLATION VALUES WITH BEST ASSAY/PREDICTOR COMBINATION (Fig. 4E) ###


message("Inferred global DNA methylation levels (Figure 4.e)...")

# find the best-performing assay in prediction of global target values via random forest:
curPredictor <- "Random Forest"
curTarget <- "Global target"
tmp <- as.data.frame(predictorResults,row.names=NA)
tmp <- tmp[tmp$metricType=="rmse" & tmp$predictor==curPredictor & tmp$targetType==curTarget,]
tmp <- tmp[order(as.numeric(as.character(tmp$value))),]
curDataset <- tmp[1,"dataset"] 

# prepare data for plotting:
ggData <- data.frame(allPredictions)
ggData <- ggData[ggData$predictor==curPredictor & ggData$targetType==curTarget & ggData$dataset==curDataset,c("dataset",sampleNames)] 
ggData[,sampleNames] <- as.numeric(ggData[,sampleNames])
ggData[,sampleNames] <- ggData[,sampleNames]/max(ggData[,sampleNames],na.rm=T) * 100
ggData <- melt(ggData)
ggData$color <- plotColLookup$sampleName[ggData$variable] 

# plot predicted DNA methylation levels:
d <- ggplot(ggData, aes(x=variable,y=value,fill=color)) + geom_bar(width=0.5,colour="#333333",size=0.25,stat="identity")  + geom_hline(yintercept=0) + scale_fill_identity() + ylab(DNA_METH_LABEL) + xlab(NULL) + defaultPlotTheme(flipX=TRUE)  + theme(axis.line.x=element_blank())
svgPlotGG(d, paste0("values_global_estimates_",curPredictor,"_",gsub(":","",gsub("/","-",curTarget)),"_",curDataset), 11, 5)



save.image(paste("methBench","afterFig4",format(Sys.time(), "%Y%m%d_%H%M"),"session.Rdata",sep="_"))
