#### ANNOTATIONS ####

# load annotations from metadata files:
datasetTable <- read.table("data/assay_table.txt",header=TRUE,sep="\t",comment.char="",quote="",na.strings="NA")
sampleTable <- read.table("data/sample_table.txt",header=TRUE,sep="\t",comment.char="",quote="",na.strings="NA")
regionTable <- read.table("data/region_table.txt",header=TRUE,sep="\t",comment.char="",quote="",na.strings="NA")
rownames(datasetTable) <- datasetTable$datasetName
rownames(sampleTable) <- sampleTable$sampleName

# prepare additional lookup tables for convenience
sampleNameLookup <- structure(c(as.character(sampleTable$sampleName),as.character(sampleTable$sampleName)), names=c(as.character(sampleTable$identifier),as.character(sampleTable$shortIdentifier)))
regionNameLookup <- structure(regionTable$regionName, names=regionTable$locus)

datasetsByType <- list()
for(curAssayType in unique(datasetTable$assayType)) {
	datasetsByType[[curAssayType]] <- rownames(datasetTable)[datasetTable$assayType==curAssayType]
}
datasetsByType[["absoluteAndRelative"]] <- union(datasetsByType[["absolute"]],datasetsByType[["relative"]])
datasetsByType[["all"]] <- rownames(datasetTable)

sampleNamesByType <- list()
for(sampleType in unique(sampleTable$sampleType)) {
	sampleNamesByType[[sampleType]] <- sampleTable[sampleTable$sampleType==sampleType,"sampleName"]
}

# select assays and regions that are to be included in the analysis:
datasetNames <- datasetsByType[["absoluteAndRelative"]]
sampleNames <- sampleTable$sampleName
regionNames <- regionTable$regionName[1:48]

# remove the region containing a SNP in a subset of samples:
regionNames <- setdiff(regionNames,REMOVE_REGIONS) 
coreRegionNames <- intersect(regionNames,paste("region_",formatC(1:16, width = 2, format = "d", flag = "0"),sep=""))

datasetNamesInOrder <- c(datasetsByType[["absolute"]],datasetsByType[["relative"]])








#### RAW DATA ####

# load locus-specific DNA methylation data:
meth <- list()
for (curDataset in datasetNames) {
	filename <- datasetTable[curDataset,"filename"]
	cat(paste("Loading",filename,"...\n"))

	# load external data file
	temp <- utils::read.table(paste0("data/",filename),header=TRUE,sep="\t",comment.char="",quote="",na.strings=MISSING_VALUES)

	# sanitise names
	names(temp) <- gsub("^X","P",names(temp))
	names(temp) <- gsub("KG1_K[0-9x]_","KG1_K_",names(temp),perl=TRUE)
	names(temp) <- gsub("KG1_T[0-9x]_","KG1_T_",names(temp),perl=TRUE)
	names(temp) <- gsub("KG1a_K[0-9x]_","KG1a_K_",names(temp),perl=TRUE)
	names(temp) <- gsub("KG1a_T[0-9x]_","KG1a_T_",names(temp),perl=TRUE)

	# pick out columns of interest
	methCols <- names(temp)[grep("*_meth",names(temp))]
	confCols <- names(temp)[grep("*_conf",names(temp))]
	otherCols <- !(names(temp)%in%methCols | names(temp)%in%confCols)
	names(temp)[otherCols] <- toCamelCase(names(temp)[otherCols])


	## dataset-specific fixes ##

	if (curDataset == "EnrichmentBS_1") {
		# take average of methylation measures on forward/reverse strand
		splitData <- split(temp[,methCols], temp$locusIdentifier)	
		combinedData <- lapply(splitData, function(X) { meanOrConcatenate(X) })		
		fwRvAvg <- do.call(rbind, lapply(combinedData, data.frame))

		# take average of averages of confidence measurements
		temp[,confCols] <- apply(temp[,confCols],c(1,2),function(x) { ifelse(is.na(x),NA,mean(as.numeric(unlist(strsplit(x,","))))) })
		splitData <- split(temp[,confCols], temp$locusIdentifier)	
		combinedData <- lapply(splitData, function(X) { meanOrConcatenate(X) })
		avgOfAvgs <- do.call(rbind, lapply(combinedData, data.frame))

		# add both together
		temp <- data.frame(locusIdentifier=row.names(fwRvAvg),fwRvAvg,avgOfAvgs)
	}
	else if (curDataset == "MS_MCA") {
		# recode qualitative values as approximate quantitative values
		temp[,methCols] <- recode(temp[,methCols],"'0.00'=0;'<0.25'=0.125;'<0,25'=0.125;'0.25-0.75'=0.5;'0,25-0,75'=0.5;'0.75-1'=0.875;'0,75-1'=0.875")
		temp[,methCols] <- apply(temp[,methCols],2,as.numeric)
	}
	else if (curDataset == "MS_HRM") {
		# recode qualitative values as approximate quantitative values
		temp[,methCols] <- apply(temp[,methCols],2,function(X) { as.numeric(recode(X,"'0'=0;'0(H)'=0;'0-0.1'=0.05;'<0.1'=0.05;'<0.1(H)'=0.05;'0.1'=0.1;'0.1(H)'=0.1;'H(0.1)'=0.1;'>0.1'=0.5;'0.1-1'=0.5;'0.1-1(H)'=0.5;'0-1'=0.5;'>0.1(H)'=0.5;'<1'=0.5;'<1(H)'=0.5;'1'=1;'H(1)'=1;'1(H)'=1;'>1'=5;'>1(H)'=5;'0-10'=5;'1-10'=5;'1-10(H)'=5;'<10'=5;'<10(H)'=5;'10'=10;'10(H)'=10;'>10'=50;'10-100'=50;'<100'=50;'100'=100;'H'=NA")) })
	}	 
	else if (curDataset == "MethyLight") {
		# replace for odd values above 100
		temp[,methCols] <- apply(temp[,methCols],2,function(x) { ifelse(x>100,100,x) })
	}
	else if (curDataset %in% c("qMSP_preamp","qMSP_standard")) {
		 # scale PCR values into an interval from zero to one
		minVal <- quantile(temp[,methCols],p=0.05,na.rm=T)
		maxVal <- quantile(temp[,methCols],p=0.95,na.rm=T)
		temp[,methCols] <- apply(temp[,methCols],2,function(x) { ifelse(x>maxVal,maxVal,x) })
		temp[,methCols] <- apply(temp[,methCols],2,function(x) { ifelse(x<minVal,0,x) })
		temp[,methCols] <- (temp[,methCols]-0)/(maxVal-0)
	}

	cat(paste("\tmin=",min(temp[,methCols],na.rm=T),", max=",max(temp[,methCols],na.rm=T),"\n",sep=""))

	# rescale [0,1] to [0,100]
	if (max(temp[,methCols],na.rm=T) <= 1) {
		cat(paste("\tAdjusting scale for",filename,"\n"))
		temp[,methCols] <- temp[,methCols]*100
	}

	if(length(meth)==0) {
		cols <- names(temp)
	} 
	else {
		newCols <- setdiff(names(temp),cols)
		if(length(newCols)>0) { cat(paste("\tUnexpected column names:",paste(newCols,collapse=","),"\n")) }
	}

	# reformat DNA methylation table
	if (length(setdiff(temp$locusIdentifier,names(regionNameLookup)))>0) {
		cat("\tIncorrect region names: ",paste(setdiff(temp$locusIdentifier,names(regionNameLookup)),collapse=","),"\n") 
	}
	names(temp) <- ifelse(names(temp)%in%names(sampleNameLookup),sampleNameLookup[names(temp)],names(temp))
	curMeth <- data.frame(regionName=regionNameLookup[temp$locusIdentifier],temp[,sampleNames])
	curMeth <- curMeth[curMeth$regionName%in%regionNames,]
	meth[[curDataset]] <- curMeth


	if(length(setdiff(names(meth[[1]]),names(temp)))>0) {
		cat(paste("\tMissing columns: ",paste(setdiff(names(meth[[1]]),names(temp)),collapse=","),"\n"))
	}
}



# prepare data for LSA plots (locus x sample x assay):
lsaData <- numeric(0)
allMeth <- data.frame(expand.grid(sampleNames,regionNames))
colnames(allMeth) <- c("sample","region")
for(curDataset in datasetNames) {
	temp <- createScatterplotTable(meth[[curDataset]])
	temp$datasetName <- curDataset
	lsaData <- rbind(lsaData,temp)
	allMeth[,curDataset] <- NA

	for(curSample in sampleNames) {
		for(curRegion in regionNames) {
			v <- temp$methValue[temp$regionName==curRegion & temp$sampleName==curSample]
			if(!is.null(v) & length(v)>0) {
				allMeth[allMeth$sample==curSample & allMeth$region==curRegion,curDataset] <- v
			}
		}
	}
}
lsaData <- lsaData[!is.na(lsaData$methValue),]
n <- nrow(lsaData)
lsaData <- data.frame(lsaData,min=rep(NA,n),max=rep(NA,n),median=rep(NA,n),mean=rep(NA,n),lower=rep(NA,n),upper=rep(NA,n))



# determine the "true" DNA methylation level for each sample/region combination as the smallest interval
# spanned by at least three differnt assay types ("consensus corridor"):
consensus <- list()
consensusFreq <- numeric(0)
for (val in c("min","max","median","mean","lower","upper")) {
	consensus[[val]] <- matrix(rep(NA,length(regionNames)*length(sampleNames)),nrow=length(regionNames),ncol=length(sampleNames),dimnames=list(regionNames,sampleNames))
}
cat("Calculating consensus corridors...\n")
for (curRegion in regionNames) {
	cat("\t* ",curRegion,"\n")
	curRegionData <- lsaData[lsaData$regionName==curRegion & lsaData$datasetName %in% datasetsByType$absolute,]
	curRegionData$assayGroup <- datasetTable[curRegionData$datasetName,"assayGroup"]
	for (curSample in sampleNames) {
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
			if (length(unique(selData$assayGroup))>=CORRIDOR_MIN_NUM_ASSAY_TYPES) break
		}

		# if we managed to find a valid interval, that's the consensus corridor:
		if (length(unique(selData$assayGroup))>=CORRIDOR_MIN_NUM_ASSAY_TYPES) {			
			consensus[["min"]][curRegion,curSample] <- min(selData$methValue,na.rm=T)
			consensus[["max"]][curRegion,curSample] <- max(selData$methValue,na.rm=T)
			consensus[["median"]][curRegion,curSample] <- median(selData$methValue,na.rm=T)
			consensus[["mean"]][curRegion,curSample] <- mean(selData$methValue,na.rm=T)			
			consensus[["lower"]][curRegion,curSample] <- max(0,min(selData$methValue,na.rm=T)-CORRIDOR_EXTENSION/2)
			consensus[["upper"]][curRegion,curSample] <- min(100,max(selData$methValue,na.rm=T)+CORRIDOR_EXTENSION/2)
			for (curCol in names(consensus)) {
				lsaData[lsaData$regionName==curRegion & lsaData$sampleName==curSample,curCol] <- consensus[[curCol]][curRegion,curSample]
			}

			consensusFreq <- rbind(consensusFreq, data.frame(
				"regionName"=curRegion,
				"sampleName"=curSample,
				"datasetName"=selData$datasetName,
				"assayGroup"=selData$assayGroup
			))
		}
		# otherwise there's no consensus value for the current region!
		# (the consensus for these will be NA, and those regions will be taken out from consensus-based analyses)
		else cat("\t\t-> No consensus found for", curRegion, " in ", curSample,": Not enough data points available.\n")
	}		
}
# calculate the deviation of each dataset's measurements from the borders of the consensus corridor:
lsaData$deviationMedian <- lsaData$methValue - lsaData$median
lsaData$deviationCorridor <- ifelse(lsaData$methValue>=lsaData$lower & lsaData$methValue<=lsaData$upper,0,ifelse(lsaData$methValue<lsaData$lower, lsaData$methValue-lsaData$lower, lsaData$methValue-lsaData$upper))



# prepare correlation tables and pair-wise dataset-vs-dataset tables:
cat("Calculating correlations and statistics...\n")
corTable <- matrix(rep(NA,length(datasetNames)^2),nrow=length(datasetNames),ncol=length(datasetNames),dimnames=list(datasetNames,datasetNames))
assayVsAssayData <- data.frame()
for(s1 in datasetNames) {
	cat("\t*",s1,"\n")
	d1 <- createScatterplotTable(meth[[s1]])
	ds1 <- melt(meth[[s1]],id.vars="regionName")
	rownames(ds1) <- paste(ds1$regionName,ds1$variable)
	for(s2 in datasetNames) {
		if (s1 != s2) {
			# prepare data
			d2 <- createScatterplotTable(meth[[s2]])
			dm <- merge(d1,d2,by=setdiff(names(d1),"methValue"))
			names(dm)[which(names(dm)=="methValue.x")] <- s1
			names(dm)[which(names(dm)=="methValue.y")] <- s2
			row.names(dm) <- dm$labelCombined
			dm <- dm[!is.na(dm[,s1])&!is.na(dm[,s2]),]
			temp <- dm[,c(s1,s2)]

			# collect statistics
			corTable[s1,s2] <- cor(temp,use="pairwise",method="pearson")[1,2]			
		} 
		else {
			corTable[s1,s2] <- 1
		}

		
		ds2 <- melt(meth[[s2]],id.vars="regionName")
		rownames(ds2) <- paste(ds2$regionName,ds2$variable)
		ds12 <- merge(ds1,ds2,by="row.names",all=F)
		ds12 <- cbind(rep(s1,nrow(ds12)),rep(s2,nrow(ds12)),ds12[,c('regionName.x','variable.x','value.x','value.y')])
		assayVsAssayData <- rbind(assayVsAssayData,ds12)
	}
}
colnames(assayVsAssayData) <- c('assay1','assay2','region','sample','meth1','meth2')
assayVsAssayData <- assayVsAssayData[complete.cases(assayVsAssayData),]


# prepare data for regression methods and classifiers:
selDatasets <- datasetsByType$absoluteAndRelative
selRegions <- coreRegionNames

mlData <- as.data.frame(matrix(NA,nrow=0,ncol=2+length(selRegions),dimnames=list(c(),c("sample","dataset",selRegions)))) #numeric(0)
for(curSample in sampleNames) {
	temp <- t(allMeth[allMeth$sample==curSample,selDatasets])
	colnames(temp) <- allMeth[allMeth$sample==curSample,"region"]
	temp <- data.frame("dataset"=rownames(temp), temp)
	mlData <- rbind(
		mlData,
		cbind("sample"=curSample,temp)	
	)
}
mlData <- mlData[,c("dataset","sample",as.character(selRegions))]
# number of valid rows per dataset:
# sapply(split(mlData,mlData$dataset),function(x) { sum(rowSums(!is.na(x))>0) })
mlData <- mlData[rowSums(!is.na(mlData[,3:ncol(mlData)]))>0,] # discard rows without any valid measurements



selRegions <- regionNames


# prepare adjusted titration data:
# (in order to account for different basal DNA methylation offsets in the titration series data,
#  we fit a linear model per region, sample, and dataset and subtract the intersect from all values.
#  Since some regions decrease in DNA methylation level (rather than increase) with increasing 
#  titration level, we also standardize the direction of change to simplify analysis and interpretation
#  by inverting the regions, for which the consensus methylation level indicates that they are anti-
#  correlated in their direction of change)
cat("Calculating offset-adjusted titration data...\n")

selRegions <- regionNames
selSamples <- sampleNames[grep("Titration",sampleNames)]
selDatasets <- datasetNames
selTitrations <- c("Titration1","Titration2")

# first determine the expected direction of change for each region/titration combination:
titrationDirections <- list()
for (titrationType in selTitrations) {
	cat("Checking expected titration directions for",titrationType,"...\n")
	nSkipped <- 0
	nSwapped <- 0
	
	titrationSamples <- selSamples[grep(titrationType,selSamples)]
	titrationPercent <- as.numeric(gsub("^001","0.01",gsub("^01","0.1",sapply(strsplit(titrationSamples,"_"),function(X) { X[3] }))))

	temp <- rep(0,length(selRegions))
	names(temp) <- selRegions
	for(curRegion in selRegions) {
		# since "true" DNA methylation levels are unknown, we take the consensus corridor medians in the extrema of the
		# titration series as a proxy to the true values and calculate the difference between these values as the "true"
		# direction of change:
		consensusMedians <- consensus[["median"]][curRegion,titrationSamples] 
		r <- cor(consensusMedians,titrationPercent,use="pairwise")
		delta <- consensusMedians[which(titrationPercent==max(titrationPercent))] - consensusMedians[which(titrationPercent==min(titrationPercent))]
		# by correlation:
		#dir <- ifelse(is.na(r) | abs(r)<0.5,0,sign(r))
		# by difference in extremes:
		dir <- ifelse(is.na(delta) | abs(delta)<5,0,sign(delta))
		cat("\t", titrationType, "-", curRegion, ":\tmeth = [", round(consensusMedians[order(titrationPercent)],1), "]\tr =", round(r,2), ", delta =", round(delta,1), "\n")

		temp[curRegion] <- dir
		if(dir==0) {
			nSkipped <- nSkipped + 1
			cat("\t\t-> Will skip",curRegion,"because there is no notable change in DNA methylation\n")
		}
		else if(dir==-1) {
			nSwapped <- nSwapped + 1
			cat("\t\t-> Will swap direction for",curRegion,"\n")
		}
					
	}
	cat("nSkipped =",nSkipped,", nSwapped =", nSwapped, "\n")

	titrationDirections[[titrationType]] <- temp
}


selData <- lsaData[lsaData$regionName%in%selRegions & lsaData$sampleName%in%selSamples,]
performanceTableTitrations <- numeric(0)

correctedTitrationData <- numeric(0)
for (titrationType in selTitrations) {
	for (curDataset in selDatasets) {
		cat("\t*",titrationType,curDataset,"\n")

		# pick out the current data:
		curSelData <- selData[selData$datasetName==curDataset & !is.na(selData$methValue) & !is.na(selData$median),]
		curSelData$titrationPercent <- as.numeric(gsub("^001","0.01",gsub("^01","0.1",sapply(strsplit(curSelData$sampleName,"_"),function(X) { X[3] }))))
		titration <- curSelData[grep(titrationType,curSelData$sampleName),]
		titrationValues <- sort(unique(titration$titrationPercent))
		if (nrow(titration)==0 | length(titrationValues)!=6) {
			cat("\t\tMissing data for",curDataset,titrationType,"\n")
			next
		}

		# for each region separately:
		for (regionName in selRegions) {
			titrationDir <- titrationDirections[[titrationType]][regionName]

			if(!is.na(titrationDir) & titrationDir!=0) { # skip those where the direction could not be determined, because consensus values were missing

				# fit a simple linear model to the data for the current assay/region/titration-series:
				curModel <- "methValue ~ x + 1"
				curData <- titration[titration$regionName == regionName,]
				if(titrationDir == -1) {
					# artificially swap direction:
					curData$methValue <- 100-curData$methValue
				}
				curData$x <- curData$titrationPercent
				if (nrow(curData)!=6) next
				fit <- lm(as.formula(curModel),data=curData)
				curData$prediction <- fitted(fit)
				f <- suppressWarnings(summary(fit))$fstatistic

				# save the corrected titraton data:
				correctedTitrationData <- rbind(correctedTitrationData, data.frame(
					"regionName" = regionName,
					"datasetName" = curDataset,
					"titrationType" = titrationType,
					"titrationPercent"=curData$titrationPercent,
					"meth"=curData$methValue,
					"methAdj"=curData$methValue-fit$coefficients["(Intercept)"]
				))
			

				# calculate performance metrics for how well the data fits the linear model
				# (and therefore how well the assay captures the linearity of the titration series)
				# and also the Pearson correlation between the observed and expeced values:
				vals <- list()
				vals[["datasetName"]] <- curDataset
				vals[["titrationType"]] <- titrationType
				vals[["regionName"]] <- regionName
				vals[["LinModel"]] <- curModel
				corTarget <- ifelse(rep(titrationType=="Titration2",length(curData$titrationPercent)),log10(curData$titrationPercent+1),curData$titrationPercent) # the "expected" target methylation level is proportional to the titration percentage!
				vals[["Pearson's r"]] <- suppressWarnings(cor(curData$methValue,corTarget))
				vals[["Adjusted R^2"]] <- suppressWarnings(summary(fit))$adj.r.squared
				vals[["Residual standard error"]] <- suppressWarnings(summary(fit))$sigma
				if (is.na(vals[["Residual standard error"]])) next
				performanceTableTitrations <- rbind(performanceTableTitrations,data.frame(vals))
			}
		}
	}
}



# load global DNA methylation data:
cat("Loading global DNA methylation data...\n")
globalMeth <- utils::read.table("data/results_global_methylation_all_assays.txt",header=TRUE,sep="\t",comment.char="",quote="",na.strings="NA",row.names=1)
# select the columns containing the (high-quality) data:
globalMethLevelsRaw <- globalMeth[,c("GlobalMeth_HPLC_MS_mean","GlobalMeth_Immunoquant_mean","GlobalMeth_Pyroseq_AluYb8_last3","GlobalMeth_Pyroseq_D4Z4_first3","GlobalMeth_Pyroseq_LINE1_first2","GlobalMeth_Pyroseq_NBL2_first4")]
# scale the values to [0,100]:
globalMethLevelsRaw <- apply(globalMethLevelsRaw,2,scaleByMax)
rownames(globalMethLevelsRaw) <- sampleNameLookup[as.character(rownames(globalMethLevelsRaw))]
colnames(globalMethLevelsRaw) <- datasetTable[as.character(datasetsByType$global),"prettyLabel"]
globalMethLevelsRaw <- globalMethLevelsRaw[sampleNames,]







#### STANDARDIZED COLORS ####

# look-up tables for colours and shapes (used in plots):
plotColLookup <- list(
	"sampleName" = structure(
		c(
			rev(rep(brewer.pal(8,"Paired")[5:6],length(sampleNamesByType$TumorNormal)/2)),
			rev(rep(brewer.pal(8,"Paired")[3:4],length(sampleNamesByType$DrugControl)/2)),
			c(
				rev(colorRampPalette(brewer.pal(8,"Paired")[1:2])(length(sampleNamesByType$Titration1))),
				rev(colorRampPalette(brewer.pal(8,"Paired")[1:2])(length(sampleNamesByType$Titration2)))
			),
			rev(rep(brewer.pal(8,"Paired")[7:8],length(sampleNamesByType$FrozenFFPE)/2))
		),
		names=as.character(c(sampleNamesByType$TumorNormal,sampleNamesByType$DrugControl,sampleNamesByType$Titration1,sampleNamesByType$Titration2,sampleNamesByType$FrozenFFPE))
	),
	"assayGroup" = structure(brewer.pal(length(unique(datasetTable[datasetNames,"assayGroup"])),"Set1"),names=as.character(unique(datasetTable[datasetNames,"assayGroup"])))
)
plotColLookup$datasetName <- structure(plotColLookup$assayGroup[datasetTable[datasetNames,"assayGroup"]],names=datasetNames)
plotColLookup$prettyLabel <- structure(plotColLookup$assayGroup[datasetTable[datasetNames,"assayGroup"]],names=datasetTable[datasetNames,"prettyLabel"])
plotColLookup$sampleType <- structure(sapply(unique(sampleTable$sampleType),function(y) { plotColLookup$sampleName[sampleTable[sampleTable$sampleType==y,"sampleName"][1]] }),names=unique(sampleTable$sampleType))


plotPchLookup <- list()
for(luType in names(plotColLookup)) {
	plotPchLookup[[luType]] <- as.numeric(as.factor(plotColLookup[[luType]]))
}








# tidy up:
suppressWarnings(rm(cutoffs,curRegionSampleData,curRegionData,curMeth,curModel,curCol,d1,d2,delta,dir,dm,ds1,ds2,ds12,forceRange,fit,fwRvAvg,i,lib,luType,n,nSkipped,nSwapped,otherCols,methCols,confCols,r,s1,s2,splitData,temp,temp1,temp2,titrationSamples,titrationDirections,titrationDir,titration,titrationValues,v,val,vals))
gc(verbose=F)

