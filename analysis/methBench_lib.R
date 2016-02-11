# Transform a label to camelCase by replacing underscores and dots
# and capitalizing the beginnings of words:
toCamelCase <- function(x){
	# inspired by http://stackoverflow.com/questions/11672050/how-to-convert-not-camel-case-to-camelcase-in-r
	capit <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2, nchar(x)))
	x <- sapply(strsplit(x, "\\."), function(x) paste(capit(x), collapse=""))
	x <- sapply(strsplit(x, "\\_"), function(x) paste(capit(x), collapse=""))
	x <- paste0(tolower(substring(x, 1, 1)), substring(x, 2, nchar(x)))
	return(x)
}

read.table <- function(...) {
	tbl <- utils::read.table(...)
	colnames(tbl) <- toCamelCase(colnames(tbl))
	return(tbl)
}

# Return the mean of all numeric columns and the concatenation of all 
# non-numeric columns:
meanOrConcatenate <- function(rows, digits=NA) {
	result <- rows[1,]
	numCols <- sapply(rows, is.numeric)
	if (is.na(digits)) {
		result[numCols] <- apply(rows[,numCols,drop=F],2,mean) 
	}
	else {
		result[numCols] <- signif(apply(rows[,numCols,drop=F],2,mean),digits)  
	}
	result[!numCols] <- sapply(rows[,!numCols,drop=F], function(X) { 
					if (length(unique(X)) == 1) return(X[1]) 
					else return(paste(X,collapse=", "))
				   })
	return(result)
}


# Reformat a ( broad)data value matrix into a long vector with 
# informative annotations:
createScatterplotTable <- function(curMeth,selRows=regionNames,selCols=sampleNames) {
	t1 <- curMeth[curMeth$regionName%in%selRows,c("regionName",selCols)]
	d1 <- expand.grid(t1$regionName,setdiff(names(t1),"regionName"),stringsAsFactors=F)
	names(d1) <- c("regionName","sampleName")
	d1 <- data.frame(d1,combinedName=apply(d1,1,paste,collapse="|"),methValue=unlist(t1[,2:ncol(t1)]))
	d1 <- data.frame(d1,labelBySample=gsub("_meth|_|","",d1$sampleName),labelByLocus=gsub("mandatory_","M",gsub("recommended_","R",gsub("optional_","O",d1$regionName))))
	d1 <- data.frame(d1,labelCombined=paste(d1$labelByLocus,d1$labelBySample,sep="|"))
	return(d1)
}

# Rescale numeric values in such a way that they fit into the 
# given range:
forceRange <- function(x,rmin,rmax) {
	if(rmin>rmax) {
		tmp <- rmin
		rmin <- rmax
		rmax <- tmp
	}
	ifelse(x>rmax,rmax,ifelse(x<rmin,rmin,x))
}

# Make a color (partially) transparent:
makeTransparent <- function(x, alpha=0.5) {
	# (see: http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color)
	if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
	alpha <- floor(255*alpha)  
	newColor <- col2rgb(col=unlist(list(x)), alpha=FALSE)
	.makeTransparent <- function(col, alpha) {
	rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
}

  newColor <- apply(newColor, 2, .makeTransparent, alpha=alpha)

  return(newColor)
}

# Recycled themes for plotting functions:
defaultPlotTheme <- function(flipX=FALSE,fontSize=10) {
	thm <- theme_bw() + theme(
				axis.line = element_line(),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				axis.title= element_text(size=fontSize,face="bold"),
				axis.text.y = element_text(hjust=1,vjust = 0.5,size = fontSize),
				axis.text.x = element_text(size = fontSize,vjust=0.5),
				strip.background = element_blank(),
				strip.text = element_text(size=fontSize*0.9,face="bold")
		) 
	if(flipX) {
		thm <- thm + theme(axis.text.x = element_text(size = fontSize, angle=90, hjust=1,vjust=0.5))
	}
	return(thm)
}
defaultTheme <- function(flipX=FALSE) {
	thm <- theme_bw() + theme(
				axis.line = element_line(),
				axis.line.y = element_line(),
				axis.line.x = element_blank(),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border = element_blank(),
				axis.title= element_text(size=8,face="bold"),
				axis.text.y = element_text(hjust=1,vjust = 0.5,size = 8),
				axis.text.x = element_text(size = 8,vjust=0.5),
				strip.background = element_blank(),
				strip.text = element_text(size=8,face="bold")
		) 
	if(flipX) {
		thm <- thm + theme(axis.text.x = element_text(size = 8, angle=90, hjust=1,vjust=0.5))
	}
	return(thm)
}


# 3-way beeswarm sample plot:
plot3WayBeeswarm <- function(d,meth="hex",extension="",dims=c("regionName","sampleName"),pch=16,bee.cex=0.65,bee.main=NA,bee.spacing=2) {
	
	for(dim1 in sort(unique(d[,dims[1]]))) {

	  	temp <- d[d[,dims[1]]==dim1,]

		sampleLevels <- rev(unique(temp[,dims[2]]))
		temp[,dims[2]] <- factor(temp[,dims[2]],levels=sampleLevels)
		
		pwcol <- unlist(plotColLookup$datasetName[temp$datasetName])

		svg(file=paste("results_analysis/plots/beeswarm_3D_",paste(dims,collapse="."),"_",extension,dim1,".svg",sep=""), width=6, height=length(sampleLevels)*0.55+1, pointsize=14)
		par(mar=c(4,7,ifelse(is.na(bee.main),3,6),1)+0.1)
		par(mgp=c(1.75,0.5,0))

	  	bees <- beeswarm(as.formula(paste("methValue","~",dims[2],sep="")), data=temp, horizontal=TRUE, cex=bee.cex, spacing=bee.spacing, cex.lab=1, pch=pch, pwcol=pwcol, bty='l', las=1,cex.axis=1,method=meth,corral="gutter",xlab=DNA_METH_LABEL,ylab="",xlim=c(-0.0001,100),main=bee.main)
		par(xpd=T) # show legend outside of the main plotting area

		assays <- unique(datasetTable[temp$datasetName,"assayGroup"])
		legend(50,length(sampleLevels)+1,assays,pch=pch,col=plotColLookup$assayGroup[assays],bty='n',horiz=FALSE,ncol=3,xjust=0.5,yjust=0)  
 
		# plot confidence interval
		corridorH <- 0.4
		for (i in 1:length(sampleLevels)) {
			dim2 <- sampleLevels[i]
			yMid <- i
			rect(xleft=consensus[["min"]][dim1,dim2], ybottom=yMid-corridorH, xright=consensus[["max"]][dim1,dim2], ytop=yMid+corridorH,col=rgb(0,0,0,0.25),border="black")
			rect(xleft=consensus[["lower"]][dim1,dim2], ybottom=yMid-0.5*corridorH, xright=consensus[["upper"]][dim1,dim2], ytop=yMid+0.5*corridorH,col=rgb(0,0,0,0.15),border="black",lwd=0.6)
	    	}

		dev.off()

	}
}

# Evaluate the results of a regression model by calculting
# the mean absolute error, root mean square error, and
# R^2 of the real and predicted data:
evaluateRegressionResults <- function(real,pred) {
	regrEvals <- regr.eval(real,pred,stats=c("mae","rmse"))

	# see: http://www.saedsayad.com/model_evaluation_r.htm
	realMean <- mean(real,na.rm=T)
	predMean <- mean(pred,na.rm=T)
	sst <- sum((real - realMean)^2, na.rm=T)
	ssr <- sum((pred - predMean)^2, na.rm=T)
	r2 <- ssr / sst

	return(c(regrEvals[1],regrEvals[2],"r2"=r2))
}
	
# Run a generalized linear model:
runGLM <- function(f, train, test, fam) {
	glmModel <- glm(f, data=train, family=fam)
	pred <- predict(glmModel, test)
	return(c("Target"=resp(f, test),"Prediction"=pred))
}
		
# Run support vector regression:
runSVM <- function(f, train, test, k) {
	require(e1071)
    	require(DMwR)

    	svmModel <- svm(f, train, type="eps-regression", scale=TRUE, kernel=k)
    	pred <- predict(svmModel, test)

	return(c("Target"=resp(f, test),"Prediction"=pred))
}		

# Run random forest regression:
runRF <- function(f, train, test) {
	require(randomForest)			
	rfModel <- randomForest(f, data=train, importance=T)
	pred <- predict(rfModel, test, type="response")
	return(c("Target"=resp(f, test),"Prediction"=pred))
}

# Scale a set of number to 0...max:
scaleByMax <- function(x,rangeUpperLimit=100) {
	x/max(x,na.rm=T)*rangeUpperLimit
}

# Impute missing values by copying the values from the nearest neighbor.
# The nearest neighbour is the most highly correlated other column.
# One particular comparison (r,c) can be disregarded for the imputation.
imputeValsFromNN <- function(x,cors=cor(x,use="pairwise",method="pearson"),r,c) {
	vals <- x[,r]
	exclIndices <- c(r,c)
	missingVals <- is.na(vals)
	nIter <- 0
	while(sum(missingVals) != 0 & nIter < ncol(cors)) {
		imputeFrom <-  names(which(cors[r,]==max(cors[r,setdiff(colnames(cors),exclIndices)])))[1]
		vals[missingVals] <- x[missingVals,imputeFrom]
		missingVals <- is.na(vals)
		exclIndices <- c(exclIndices,imputeFrom)
		nIter <- nIter + 1
	}
	if(nIter==ncol(cors)) {
		message("Warning: Imputation was impossible for: ", r, "x", c," (tried:",paste(exclIndices,collapse="/"),")\n")
	}
	return(vals)
}

# Calculate the correlation table between a numeric matrix. In doing so,
# impute missing values by copying them from the most similar (correlated)
# column that is not being looked at in the calculation of the current
# correlation coefficient:
nnImputedCor <- function(x,method="pearson") {
	tmpCorTable <- cor(x,use="pairwise",method=method)
	corTable <- matrix(NA,ncol=ncol(x),nrow=ncol(x))
	colnames(corTable) <- rownames(corTable) <- colnames(x)

	for(r in colnames(x)) {
		for(c in colnames(x)) {		
			ds1 <- imputeValsFromNN(x,tmpCorTable,r,c)
			ds2 <- imputeValsFromNN(x,tmpCorTable,c,r)
			corTable[r,c] <- cor(ds1,ds2,method=method)			
			#message("* cor:",r,c,"=",corTable[r,c],"\n")
		}
	}

	return(corTable)
}

# Make a set of target values for the prediction by regression models
# by taking assorted real data and removing outliers:
makeGlobalReference <- function(x,useData) {
	refData <- x[,useData]
	# remove titration outliers by imputing the value as the mean of the surrounding titration levels:
	for(tNum in 1:2) {
		tType <- paste("Titration",tNum,sep="")
		selSamples <- sort(sampleNames[grepl(tType,sampleNames)])
		for(sampleI in 2:(length(selSamples)-1)) {
			tmp <- x[selSamples[c(sampleI-1,sampleI,sampleI+1)],]
			for(curDataset in useData) {
				if(tmp[2,curDataset]>tmp[1,curDataset]) {
					#message("replace ", curDataset, ": ", rownames(tmp),"\n")
					refData[selSamples[sampleI],curDataset] <- mean(tmp[1,curDataset],tmp[3,curDataset])
				}
			}
		}
		# last one in series:
		sampleI <- length(selSamples)
		tmp <- x[selSamples[c(sampleI-1,sampleI)],]
		for(curDataset in useData) {
			if(tmp[2,curDataset]>tmp[1,curDataset]) {
				#message("replace last one in series ", curDataset, ": ", rownames(tmp),"\n")
				refData[selSamples[sampleI],curDataset] <- tmp[1,curDataset]
			}
		}
	}
	apply(refData,1,mean)
}

# Calculate feature importance (for each region) by using 
# an F-score:
importanceByFScore <- function(data) {
	# see Chang YW & Lin CJ. Feature Ranking Using Linear SVM. JMLR: Workshop and Conference Proceedings 3: CI2008 workshop on causality (2008).
	r <- intersect(selRegions,colnames(data))
	fg <- data[,"Target"]==unique(data[,"Target"])[1]
	bg <- data[,"Target"]==unique(data[,"Target"])[2]
	x <- data[,r]
	xPlus <- x[fg,]
	xMinus <- x[bg,]
	
	sqDiffMean <- function(x){ sum(x - mean(x))^2 }

	f <- ((colMeans(xPlus) - colMeans(x))^2 + (colMeans(xMinus) - colMeans(x))^2)
	f <- f / ( ((1 / (nrow(xPlus)-1)) * apply(xPlus,2,sqDiffMean)) + ((1 / (nrow(xMinus)-1)) * apply(xMinus,2,sqDiffMean)) ) 
	f
}

# Generate a plot visualizing the change in ROC area under curve 
# upon changing some sort of parameter in the model training, e.g.
# the noise level:
plotAUCChange <- function(d,by,categoryTypesName,categoryTypes=unique(d[,categoryTypesName]),selectedDatasets=unique(d[,"datasetName"]),exclude=c(),asc=F) {
	ggData <- as.data.frame(d)
	ggData <- ggData[!is.na(ggData[,by]),]
	ggData[,by] <- factor(as.character(ggData[,by]),levels=as.character(sort(unique(ggData[,by]),decreasing=!asc)))
	par(mfrow=c(1,2))
	par(mar=c(1,4,9,1)+0.1)
	par(mgp=c(2,0.65,0))
	orderedBy <- as.character(levels(ggData[,by]))
	orderedByTidy <- setdiff(orderedBy,as.character(exclude))
	n <- length(orderedBy)
	syms <- c(17,19,15)
	if(length(orderedByTidy)!=3) {
		syms <- 1:length(orderedByTidy)
	}
	pchs <- structure(syms,names=orderedByTidy)

	arrowMargin <- 0.015

	for(categoryType in categoryTypes) {
		x <- 0.5:(length(selectedDatasets)+0.5)
		yRange <- c(0.3,1) # c(0.75*min(as.numeric(ggData$AUC)),1)
		plot(x,rep(1,length(x)),type="l",bty="n",ylim=yRange,xlim=range(x),xaxt="n",xlab=NA,ylab="ROC area under curve",col="white")
		axis(3,1:length(selectedDatasets),labels=datasetTable[selectedDatasets,"prettyLabel"],las=3,tick=F)
		abline(h=1,col="black")
		abline(h=0.5,col="darkgray",lty=2)
		for(curDataset in selectedDatasets) { 
			dsIndex <- which(curDataset==selectedDatasets)
			tmp <- ggData$dataset==curDataset & ggData[,categoryTypesName]==categoryType
			tmp <- structure(ggData[tmp,"AUC"],names=as.character(ggData[tmp,by]))

			plotCol <- plotColLookup$datasetName[curDataset]
			if(curDataset %in% datasetsByType$relative) {
				plotCol <- scales::muted(plotCol, l=85, c=90)
			}
			tmp <- as.numeric(tmp[orderedByTidy])

			points(rep(dsIndex,length(orderedByTidy)),tmp,col=plotCol,pch=pchs[orderedByTidy],cex=1.15)
			arrFrom <- tmp[1:(length(tmp)-1)]
			arrTo <-  tmp[2:length(tmp)]
			gtMinLen <- abs(arrFrom-arrTo)>(3*arrowMargin)
			if(sum(gtMinLen)>0) {
				suppressWarnings(arrows(dsIndex, (pmax(arrFrom,arrTo) - arrowMargin)[gtMinLen], dsIndex, (pmin(arrFrom,arrTo) + arrowMargin)[gtMinLen],col=plotCol,length=0,lwd=0.5))
			}
		}
		legend("bottom",orderedByTidy,pch=pchs[orderedByTidy], bty="n", horiz=TRUE, title=categoryType)
	}
}

# Add a random error to the given measurements by replacing a
# percentage of real values with random values drawn from the same
# distribution:
addRandomError <- function(x,n=0.1,universe=unlist(x)){
	tmp <- x[,2:ncol(x)]
	rnd <- matrix(runif(ncol(tmp)*nrow(tmp),max=1,min=0),ncol=ncol(tmp),nrow=nrow(tmp))<=n 
	tmp[rnd] <- sample(universe,sum(rnd),replace=T)
	x[,2:ncol(x)] <- tmp
	return(x)
}

# Add uniform noise to the data by altering the values by adding
# or subtracting a random number:
addUniformNoise <- function(x,n=0.1,universe=unlist(x)){
	tmp <- data.matrix(x[,2:ncol(x)])
	rng <- (max(tmp,na.rm=T)-min(tmp,na.rm=T))*n
	x[,2:ncol(x)] <- tmp + matrix(runif(ncol(tmp)*nrow(tmp),max=rng,min=-rng),ncol=ncol(tmp),nrow=nrow(tmp))  
	return(x)
}

# Imput missing values from the median of values from other 
# observations, but use only the other obervations from the same
# group to do so. Otherwise what happens is this:
# If a value is missing for one sample and there's an equal number
# of samples in both groups, taking the median across all values
# biases the imputed value towards the opposite group. If there 
# are many observations, this hardly matter, however, we only have
# a very small test sample in this study! 
groupAwareImpute <- function(x) {
	# first split the data by group:
	gLabels <- unique(x[,1])
	g1 <- as.data.frame(x[x[,1]==gLabels[1],2:ncol(x)])
	g2 <- as.data.frame(x[x[,1]==gLabels[2],2:ncol(x)])

	# impute values by column median:
	g1 <- impute(g1,"median")
	g2 <- impute(g2,"median")

	# find all-NA columns:
	allNA <- colSums(!is.na(g1))==0 | colSums(!is.na(g2))==0
	if(sum(allNA)>0) {
		message("Columns needed to be removed because they were all NA in in one group: ", colnames(g1)[allNA])
		g1 <- g1[,!allNA]
		g2 <- g2[,!allNA]
	}

	rbind(data.frame("Target"=gLabels[1],g1),data.frame("Target"=gLabels[2],g2))
}

# Run a number of random trials using the given transformation to the data (rndDataModFun)
# and performing patient-stratified n-fold crossvalidation in each step:
runClassifierTrials <- function(x,rndDataModFun=identity,nTrials=CLASSIFIER_TRIALS,comparisons=cmp,testTrainCombos=combos) {
	require("e1071")

	f <- as.factor(Target) ~ .

	ps <- c()
	ts <- c()
	
	runClassifier <- function(trainD) {
		svm(f, data=trainD, probability = TRUE, scale=FALSE, kernel="linear")
	}

	x <- groupAwareImpute(x)

	nSkipped <- 0
	for(n in 1:nTrials) {
		y <- rndDataModFun(x)
		y <- y[sample(nrow(y),nrow(y)),] # need to shuffle the data so that observations for the two classes are in mixed order

		y[,2:ncol(y)] <- y[,2:ncol(y)]/100 # for SVM: put values into range [0;1]

		if(n %% 10 == 0) message(n," ")

		for(testIndex in 1:nrow(testTrainCombos)) {
			s <- as.character(testTrainCombos[testIndex,])
			trainD <- as.data.frame(y[setdiff(unlist(comparisons),s),])
			testD <- as.data.frame(y[s,])
			
			tmp <- as.matrix(rbind(trainD,testD)[,2:ncol(trainD)])
			if(sum(apply(tmp,2,sd)!=0)==0) {
				nSkipped <- nSkipped+1
				next;
			}
		
			classifierModel <- runClassifier(trainD)
			pred <- predict(classifierModel, testD, decision.values=TRUE)
			ps <- c(ps,attr(pred,"decision.values"))

			ts <- c(ts,as.numeric(testD$Target))
		}
	}
	message("")
	if(nSkipped>0) message("nSkipped =",nSkipped," (of",nrow(testTrainCombos)*nTrials,")")

	perf <- performance(prediction(ps, ts), measure = "auc", x.measure = "cutoff")
	auc <- as.numeric(deparse(as.numeric(perf@y.values)))

	perf <- performance(prediction(matrix(ps,nrow=(nTrials*nrow(testTrainCombos))-nSkipped), matrix(as.factor(ts),nrow=(nTrials*nrow(testTrainCombos))-nSkipped)), measure = "tpr", x.measure = "fpr")


	return(list(
		"auc"=auc,
		"perf"=perf
	))
}
	
# Plot selected ROC curves:
plotSelectedNoiseRocs <- function(curDataset, withLegend=FALSE) {
	tmpData <- rocPlotDataClassifiers[[curDataset]]

	first <- TRUE

	# then train classifiers with the different noise/error models...
	for(noiseType in names(NOISE_TYPES)) {	
		tmpData2 <- tmpData[[noiseType]]
		for(noiseLevelI in 1:length(setdiff(NOISE_LEVELS,1))) {
			noiseLevel <- NOISE_LEVELS[noiseLevelI]
			if(noiseLevel==0 & noiseType!=names(NOISE_TYPES)[1]) next;

			tmpData3 <- tmpData2[[as.character(noiseLevel)]]

			plot(tmpData3, col=NOISE_COLORS[[noiseType]][noiseLevelI],add=!first, lty=which(NOISE_LEVELS==noiseLevel), downsampling=500, main=datasetTable[curDataset,"prettyLabel"],lwd=2.25, avg="threshold")
			first <- F
		}
	}
	if(withLegend) {
		legend("bottomright",as.character(setdiff(NOISE_LEVELS,1)),lty=1:length(setdiff(NOISE_LEVELS,1)), bty = "n", cex = 1, title="Noise Level",lwd=2.75)
		cols <- c()
		for(noiseType in names(NOISE_TYPES)) {
			cols <- c(cols,NOISE_COLORS[[noiseType]][2])
		}
		legend("bottom",as.character(names(NOISE_TYPES)),col=cols, cex = 1, title="Noise Type", bty = "n", lwd=2.75)
	}
}
plotSelectedSelectionRocs <- function(curDataset, withLegend=FALSE) {
	tmpData <- rocPlotDataClassifiersSelection[[curDataset]]

	first <- TRUE

	# then train classifiers with the different noise/error models...
	for(regionSelectionMethodName in names(REGION_SELECTION_METHODS)) {
		tmpData2 <- tmpData[[regionSelectionMethodName]]
		for(howMany in REGION_NUMBERS) {
			tmpData3 <- tmpData2[[as.character(howMany)]]

			plot(tmpData3, col=selectionColors[which(regionSelectionMethodName==names(REGION_SELECTION_METHODS))],add=!first, lty=which(howMany==REGION_NUMBERS), main=datasetTable[curDataset,"prettyLabel"],downsampling=500, lwd=2.25, avg="threshold")	
			first <- FALSE
		}
	}
	if(withLegend) {
		legend("bottomright",as.character(REGION_NUMBERS),lty=1:length(REGION_NUMBERS), bty = "n", cex = 1, title="Number of regions",lwd=2.75)
		legend("bottom",as.character(names(REGION_SELECTION_METHODS)),col=selectionColors, cex = 1, title="Selection method", bty = "n", lwd=2.75)
	}
}


# convert DNA methylation beta- to m-values:
beta2mval <- function (betas, epsilon = 0.00001)  {
	# function extracted from RnBeads (Assenov et al., Nature Methods, 2014)
	if (!is.numeric(betas)) {
		stop("invalid value for betas")
	}
	if (!(is.numeric(epsilon) && length(epsilon) == 1 && (!is.na(epsilon)))) {
		stop("invalid value for epsilon")
	}
	if (epsilon < 0 || epsilon > 0.5) {
		stop("invalid value for epsilon; expected 0 <= epsilon <= 0.5")
	}
	betas[betas < epsilon] <- epsilon
	betas[betas > (1 - epsilon)] <- 1 - epsilon

	return(log2(betas/(1 - betas)))
}

# Plotting device functions:
svgPlot <- function(name,w,h,pointsize=12) {	
	CairoSVG(file=paste0("results_analysis/plots/",name,".svg"), width=w, height=h, pointsize=pointsize)
}
svgPlotGG <- function(p,name,w,h,units="cm") {
	suppressWarnings(ggsave(filename=paste0("results_analysis/plots/",name,".svg"), plot=p, scale=1, width=w, height=h, units=units))
}

# Pick a number of random features from the given data:
pickRandomRegions <- function(x,howMany=5) {
	tmp <- data.matrix(x[,2:ncol(x)])
	x <- data.frame(cbind("Target"=x[,"Target"],tmp[,sample(1:ncol(tmp),howMany)]))
	return(x)
}

# Sanitize sample names to make them more consistent:
sanitizeSampleNames <- function(x) {
	x <- gsub("^ ","",x)	
	x <- gsub(" $","",x)	
	x <- gsub("IVM\\-?(\\d)","IVM_\\1",x,perl=TRUE)	
	x <- gsub("^X","P",x)
	x <- gsub("KG1[\\s_](T|K).*","KG1 \\11",x,perl=TRUE)
	x <- gsub("KG1a[\\s_](T|K).*","KG1a \\11",x,perl=TRUE)
	x <- gsub("HCT115","HCT15",x,perl=TRUE)
	x <- gsub("_a$","",x,perl=TRUE)
	x <- gsub("B_(\\d)","B\\1",x,perl=TRUE)
	#x <- gsub("HCT115","HCT15",x,perl=TRUE)
	return(x)
}

# Pick a number of highly informative features from the
# given data. Information content is calculated as the
# F-score of each feature in the given data.
pickTopRegions <- function(x,howMany=1) {
	regionImportance <- names(sort(importanceByFScore(x),decreasing=T))
	as.matrix(x[,c("Target",regionImportance[1:howMany])])
}

# the R libary "car" (https://cran.r-project.org/web/packages/car/) turned out to be
# a little tricky to install on our system, so we extracted the two following functions
# we wanted to use:
recode <- function (var, recodes, as.factor.result, as.numeric.result = TRUE, levels) {
    lo <- -Inf
    hi <- Inf
    recodes <- gsub("\n|\t", " ", recodes)
    recode.list <- rev(strsplit(recodes, ";")[[1]])
    is.fac <- is.factor(var)
    if (missing(as.factor.result)) 
        as.factor.result <- is.fac
    if (is.fac) 
        var <- as.character(var)
    result <- var
    for (term in recode.list) {
        if (0 < length(grep(":", term))) {
            range <- strsplit(strsplit(term, "=")[[1]][1], ":")
            low <- try(eval(parse(text = range[[1]][1])), silent = TRUE)
            if (class(low) == "try-error") {
                stop("\n  in recode term: ", term, "\n  message: ", 
                  low)
            }
            high <- try(eval(parse(text = range[[1]][2])), silent = TRUE)
            if (class(high) == "try-error") {
                stop("\n  in recode term: ", term, "\n  message: ", 
                  high)
            }
            target <- try(eval(parse(text = strsplit(term, "=")[[1]][2])), 
                silent = TRUE)
            if (class(target) == "try-error") {
                stop("\n  in recode term: ", term, "\n  message: ", 
                  target)
            }
            result[(var >= low) & (var <= high)] <- target
        }
        else if (0 < length(grep("^else=", squeezeBlanks(term)))) {
            target <- try(eval(parse(text = strsplit(term, "=")[[1]][2])), 
                silent = TRUE)
            if (class(target) == "try-error") {
                stop("\n  in recode term: ", term, "\n  message: ", 
                  target)
            }
            result[1:length(var)] <- target
        }
        else {
            set <- try(eval(parse(text = strsplit(term, "=")[[1]][1])), 
                silent = TRUE)
            if (class(set) == "try-error") {
                stop("\n  in recode term: ", term, "\n  message: ", 
                  set)
            }
            target <- try(eval(parse(text = strsplit(term, "=")[[1]][2])), 
                silent = TRUE)
            if (class(target) == "try-error") {
                stop("\n  in recode term: ", term, "\n  message: ", 
                  target)
            }
            for (val in set) {
                if (is.na(val)) 
                  result[is.na(var)] <- target
                else result[var == val] <- target
            }
        }
    }
    if (as.factor.result) {
        result <- if (!missing(levels)) 
            factor(result, levels = levels)
        else as.factor(result)
    }
    else if (as.numeric.result && (!is.numeric(result))) {
        result.valid <- na.omit(result)
        opt <- options(warn = -1)
        result.valid <- as.numeric(result.valid)
        options(opt)
        if (!any(is.na(result.valid))) 
            result <- as.numeric(result)
    }
    result
}

# Drop white-space:
squeezeBlanks <- function (text) {
    gsub(" *", "", text)
}
