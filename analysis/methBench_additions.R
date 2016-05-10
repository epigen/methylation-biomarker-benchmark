# ADDITIONAL DATA AND ANALYSES ADDED DURING REVISION


# ensure reproducibility by fixing the random seed to an arbitrary offset:
load(latestPrepName)
set.seed(1234)
message("=== ADDITIONAL DATA AND ANALYSES ===")




##############################################################################
############ CLONAL BISULFITE SEQUENCING (Supplementary Figure 7) ############
##############################################################################

dClonalBS <- read.table("data/results_validation_ClonalBS.txt", sep="\t", header=TRUE, comment.char="", quote="", check.names=FALSE)
rownames(dClonalBS) <- dClonalBS$locusIdentifier
rownames(dClonalBS)[!is.na(regionNameLookup[dClonalBS$locusIdentifier])] <- regionNameLookup[dClonalBS$locusIdentifier][!is.na(regionNameLookup[dClonalBS$locusIdentifier])] 
dClonalBS <- dClonalBS[,grep("Meth",colnames(dClonalBS))]
colnames(dClonalBS) <- sampleNameLookup[gsub("Meth", "", colnames(dClonalBS), ignore.case=TRUE)]
selRegions <- intersect(regionNames,rownames(dClonalBS)[rowSums(!is.na(dClonalBS))>0])
selSamples <- colnames(dClonalBS)[colSums(!is.na(dClonalBS))>0]
dOther <- dcast(melt(allMeth[allMeth$region%in%selRegions & allMeth$sample%in%selSamples,]), variable + region ~ sample, value.var="value")
dOther <- lapply(split(dOther,dOther$variable), function(x) { rownames(x) <- x$region; x[selRegions,selSamples] })
dClonalBS <- dClonalBS*100.0


### REPLICATES COMPARISON (Supplementary Figure 7.a) ###

svgPlot("clonal_reps", 3, 3)
par(mar=c(4,4,1,1)+0.1)
x <- c(dClonalBS["region_07",],dClonalBS["region_08",])
y <- c(dClonalBS["mandatory_7_b",],dClonalBS["mandatory_8_b",])
plot(as.numeric(x), as.numeric(y), xlab="ClonalBS 1", ylab="ClonalBS 2", pch=16, col=rgb(0.1,0.1,0.1,0.3), xlim=c(0,100), ylim=c(0,100), bty="l")
legend("topleft", paste0("r = ", round(cor(as.numeric(x), as.numeric(y), use="pairwise", method="pearson"),2)), bty="n")
dev.off()

### CORRELATION HEATMAP (Supplementary Figure 7.b) ###

tmpSelRegions <- c("region_07","region_08")
dCor <- dOther[intersect(names(dOther),datasetsByType$absolute)]
names(dCor) <- datasetTable[names(dCor),"prettyLabel"]
dCor <- dcast(melt(c(list("ClonalBS 1"=cbind(tmpSelRegions,dClonalBS[tmpSelRegions,]), "ClonalBS 2"=cbind(tmpSelRegions,dClonalBS[c("mandatory_7_b","mandatory_8_b"),]) ),lapply(dCor, function(x) cbind(tmpSelRegions,x[tmpSelRegions,selSamples]))), id.var=c("tmpSelRegions")), L1~tmpSelRegions+variable, value.var="value")
rownames(dCor) <- dCor[,1]
dCor <- dCor[,-1]
dCor <- cor(t(dCor), use="pairwise", method="pearson")

CairoSVG("results_analysis/plots/clonal_cors_reg78", width=13, height=13, pointsize=12)
pheatmap(data.matrix(dCor), ,col=colorRampPalette(rev(brewer.pal(5,"Oranges")),bias=0.1,space="Lab")(20*5),breaks=seq(0,1,0.01),border_color="white",cellwidth=22,cellheight=22, cluster_cols=TRUE, cluster_rows=TRUE, treeheight_col=0, number_color="black", display_numbers=TRUE)
dev.off()

### DEVIATION FROM CORRIDOR (Supplementary Figure 7.c) ###

# only proceed with selected regions and samples:
dClonalBS <- dClonalBS[selRegions, selSamples]
dClonalDev <- melt(cbind("region"=rownames(dClonalBS[selRegions,selSamples]),dClonalBS[selRegions,selSamples]))
dClonalDev <- data.frame(dClonalDev, t(apply(dClonalDev,1,function(x) {
	r <- x["region"]
	s <- x["variable"]
	c("lower"=consensus$lower[r,s],"upper"=consensus$upper[r,s],"median"=consensus$median[r,s])
})))
dClonalDev$deviation <- ifelse(dClonalDev$value>=dClonalDev$lower & dClonalDev$value<=dClonalDev$upper,0,ifelse(dClonalDev$value<dClonalDev$lower, dClonalDev$value-dClonalDev$lower, dClonalDev$value-dClonalDev$upper))
dClonalDev$tn <- gsub("CRC_._","",dClonalDev$variable,perl=TRUE)
dClonalDev$crc <- substr(dClonalDev$variable,1,5)

# statistics for the supplementary table:
message("ClonalBS: mean = ", round(mean(unlist(dClonalBS), na.rm=T),3), ", sd = ", round(sd(unlist(dClonalBS), na.rm=T),3), ", min = ", round(min(unlist(dClonalBS), na.rm=T),3), ", max =", round(max(unlist(dClonalBS), na.rm=T),3))
message("ClonalBS: quantiles (0.1,0.25,0.5,0.75,0.9) = ", paste(round(quantile(unlist(dClonalBS),c(0.1,0.25,0.5,0.75,0.9), na.rm=T),3), collapse=", "))
message("ClonalBS: deviation - mean abs. = ", round(mean(abs(dClonalDev$deviation), na.rm=T),3), ", median = ", round(median(abs(dClonalDev$deviation), na.rm=T),3), ", directional = ", round(mean(dClonalDev$deviation, na.rm=T),3))
 
# absolute and directional deviation plots: 
tmp <- lsaData[ lsaData$sampleName%in%selSamples & lsaData$regionName%in%selRegions & lsaData$datasetName%in%datasetsByType$absolute,c("datasetName","deviationCorridor")]
tmp$datasetName <- datasetTable[tmp$datasetName,"prettyLabel"]
allDevData <- rbind(data.frame("datasetName"="ClonalBS","deviationCorridor"=dClonalDev[,"deviation"]),tmp)
allDevData$datasetName <- factor(allDevData$datasetName, levels=rev(c("ClonalBS",datasetTable[datasetNamesInOrder,"prettyLabel"])) )
p <- ggplot(allDevData, aes(x=datasetName, y=abs(deviationCorridor), fill=datasetName)) + geom_boxplot() + defaultPlotTheme(flipX=FALSE) + xlab(NULL) + ylab("Absolute deviation") + scale_fill_manual(values=c("ClonalBS"="gray",plotColLookup$prettyLabel), guide=FALSE) + coord_flip() + ylim(0,100)
svgPlotGG(p, "clonal_box", 9, 19)
p <- ggplot(aggregate(deviationCorridor~datasetName,allDevData,mean), aes(x=datasetName, y=deviationCorridor, fill=deviationCorridor)) + geom_bar(colour="#333333", stat="identity") + defaultPlotTheme(flipX=FALSE) + xlab(NULL) + ylab("Directional deviation") + coord_flip() + scale_fill_gradient2(low=brewer.pal(3,"BrBG")[3],high=brewer.pal(3,"BrBG")[1],guide=F,space="Lab") + geom_hline(aes(yintercept=0),colour="#333333") + ylim(-5,5)
svgPlotGG(p, "clonal_box_dir", 6.5, 19)

### SCATTER PLOT PANEL (Supplementary Figure 7.d) ###

dScatter <- dcast(rbind(data.frame(melt(cbind("region"=rownames(dClonalBS[selRegions,selSamples]),dClonalBS[selRegions,selSamples])),"L1"="ClonalBS 1"),melt(lapply(dOther, function(x) { cbind("region"=rownames(x[selRegions,selSamples]),x[selRegions,selSamples]) }))), L1~region+variable, value.var="value")
rownames(dScatter) <- dScatter[,1]
dScatter <- dScatter[,-1]
dScatterX <- numeric()
for(s2 in c("ClonalBS")) {
	for(s1 in datasetsByType$absolute) {
		tmp <- data.frame("meth1"=as.numeric(dScatter[s1,]),"meth2"=as.numeric(dScatter[s2,]), "dataset"=datasetTable[s1,"prettyLabel"],"txt"=NA)
		tmp[1,"txt"] <- paste("r =",format(cor(as.numeric(dScatter[s1,]),as.numeric(dScatter[s2,]), use="pairwise"),digits=3,nsmall=3))
		dScatterX <- rbind(tmp,dScatterX)
		
	}
}
d <- ggplot(data=dScatterX,aes(x=meth1,y=meth2)) + geom_point(shape=16,alpha=0.3) + geom_smooth(method=lm,fill='steelblue') + theme_bw() + xlab(DNA_METH_LABEL) + ylab(DNA_METH_LABEL) + xlim(0,100) + ylim(0,100) + defaultPlotTheme(fontSize=16) + theme(axis.line = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank()	) + facet_wrap(~dataset, ncol=6) + geom_hline(y_intercept=0) + geom_vline(x_intercept=0) + geom_text(aes(label=txt), x=3, y=100, ,hjust=0,vjust=1,size=4, stat = "identity", position = "identity")
svgPlotGG(d, "clonal_scatters", 12, 7, units="in")




######################################################################################
############ LOW-INPUT TITRATION SERIES (Supplementary Figures 10 and 11) ############
######################################################################################

dataDir <- "data/"
dLowInput <- numeric()
dLowInputStatus <- numeric()
dLowInputSuccess <- data.frame()

# load in and recode data:
for(f in list.files(dataDir, pattern="results_low_input")) {
	tmp <- read.table(paste0(dataDir,f), sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)
	datasetName <- gsub("results_low_input_(.+).txt", "\\1", f, perl=TRUE)

	methCols <- grep("Meth",colnames(tmp),value=TRUE)

	# transform ranges into quantitative measurements:
	if (datasetName == "MS_MCA") {
		tmp[,methCols] <- recode(as.matrix(tmp[,methCols]),"'failed'='failed';'0.00'=0;'<0.25'=0.125;'<0,25'=0.125;'0.25-0.75'=0.5;'0,25-0,75'=0.5;'0.75-1.00'=0.875;'0,75-1.00'=0.875;'0.75-100'=0.875")
	}
	else if (datasetName == "MS_HRM") {
		tmp[,methCols] <- apply(as.matrix(tmp[,methCols]),2,function(X) { recode(X,"'0'=0;'0(H)'=0;'0-0.1'=0.05;'<0.1'=0.05;'<0.1(H)'=0.05;'0.1'=0.1;'0.1(H)'=0.1;'H(0.1)'=0.1;'>0.1'=0.5;'0.1 - 1'=0.5;'0.1-1'=0.5;'0.1-1(H)'=0.5;'0-1'=0.5;'>0.1(H)'=0.5;'<1'=0.5;'<1(H)'=0.5;'1'=1;'H(1)'=1;'1(H)'=1;'>1'=5;'>1(H)'=5;'0-10'=5;'1-10'=5;'1-10(H)'=5;'<10'=5;'<10(H)'=5;'10'=10;'10(H)'=10;'>10'=50;'10-100'=50;'10-100'=50;'<100'=50;'100'=100;'H'=NA") })
	}

	# summarize status in terms of "good" or "failed":
	tmpR <- regionNameLookup[tmp$locusIdentifier]
	tmpM <- suppressWarnings(data.matrix(tmp[,methCols]))
	if(max(tmpM,na.rm=T)<=1) tmpM <- tmpM*100
	colnames(tmpM) <- tolower(gsub("Meth","",methCols))
	resultStatus <- as.matrix(tmp[,methCols])
	colnames(resultStatus) <- tolower(gsub("Meth","",methCols))
	resultStatus[is.finite(tmpM)] <- "Good"
	resultStatus[resultStatus%in%c("FAILED","FAILED ","failed","signal too low")] <- "Failed"	
	resultStatus[resultStatus%in%c("FAILED : unspecific PCR","no amplicon","failed/digital","no reads")] <- "Failed"
	resultStatus[is.na(resultStatus) | resultStatus%in%c("")] <- "N/A"
	void <- sapply(colnames(resultStatus),function(x) {
		dLowInputStatus <<- rbind( data.frame("datasetName"=datasetName,"variable"=x,data.frame(table(resultStatus[,x]))) ,dLowInputStatus)
	})

	dLowInputSuccess <- rbind(data.frame("datasetName"=datasetName,melt(resultStatus)),dLowInputSuccess)
	colnames(tmpM) <- tolower(gsub("Meth","",methCols))

	dLowInput <- rbind(melt(data.frame("datasetName"=datasetName,"regionName"=as.character(tmpR),tmpM), id.vars=c("datasetName","regionName")),dLowInput)
}
colnames(dLowInputStatus) <- c("datasetName","variable","status","freq")

# tidy and format data:
dLowInput <- dLowInput[!is.na(dLowInput$value),]
dLowInput <- dLowInput[!(dLowInput$regionName%in%c("region_11")),] # remove SNP region
dLowInput$series <- factor(ifelse(toupper(substr(dLowInput$variable,1,1))=="X","Fragmented","Intact"),levels=c("Intact","Fragmented"))
dLowInput$quantity <- as.numeric(gsub("03","0.3",gsub("ng","",substring(dLowInput$variable,2))))
dLowInput$quantity <- factor(dLowInput$quantity, levels=sort(unique(dLowInput$quantity), decreasing=TRUE))
dLowInputStatus$series <- as.character(ifelse(toupper(substr(dLowInputStatus$variable,1,1))=="X","Fragmented","Intact"),levels=c("Intact","Fragmented"))
dLowInputStatus$quantity <- as.numeric(gsub("03","0.3",gsub("ng","",substring(dLowInputStatus$variable,2))))
dLowInputStatus$quantity <- factor(dLowInputStatus$quantity, levels=sort(unique(dLowInputStatus$quantity), decreasing=TRUE))
dLowInput$prettyLabel <- factor(datasetTable[dLowInput$datasetName, "prettyLabel"], levels=datasetTable[datasetNamesInOrder[datasetNamesInOrder%in%dLowInput$datasetName], "prettyLabel"])
dLowInput$assayGroup <-  datasetTable[dLowInput$datasetName, "assayGroup"]
dLowInput$sampleName <-  sampleNameLookup[ifelse(dLowInput$series=="Fragmented", "1T", "2T")]
dLowInputStatus$prettyLabel <- factor(datasetTable[dLowInputStatus$datasetName, "prettyLabel"], levels=datasetTable[datasetNamesInOrder[datasetNamesInOrder%in%dLowInputStatus$datasetName], "prettyLabel"])
dLowInputStatus$status <- factor(dLowInputStatus$status, levels=c("Good","Failed","No quantitation"))
dLowInputStatus <- dLowInputStatus[!is.na(dLowInputStatus$status),]

# determine success rates from status values:
dLowInputSuccessRate <- merge(aggregate(freq~prettyLabel+series+quantity, dLowInputStatus[dLowInputStatus$status=="Good",], sum), aggregate(freq~prettyLabel+series+quantity, dLowInputStatus, sum), by=c("prettyLabel", "series", "quantity"), all=TRUE)
dLowInputSuccessRate$freq.x[is.na(dLowInputSuccessRate$freq.x)] <- 0
dLowInputSuccessRate$successRate <- dLowInputSuccessRate$freq.x/dLowInputSuccessRate$freq.y
dLowInputSuccessRateTotal <- merge(aggregate(freq~prettyLabel+series, dLowInputStatus[dLowInputStatus$status=="Good",], sum), aggregate(freq~prettyLabel+series, dLowInputStatus, sum), by=c("prettyLabel", "series"), all=TRUE)
dLowInputSuccessRateTotal$freq.x[is.na(dLowInputSuccessRateTotal$freq.x)] <- 0
dLowInputSuccessRateTotal$successRate <- dLowInputSuccessRateTotal$freq.x/dLowInputSuccessRateTotal$freq.y

# calculate deviation of each measurement from "own target" value, i.e. the measurement
# in the corresponding high-input experiment:
dLowInput$ownTarget <- apply(dLowInput,1,function(x) allMeth[allMeth$region==x["regionName"] & allMeth$sample==x["sampleName"], x["datasetName"]] )
dLowInput$ownTargetDeviation <- dLowInput$value-dLowInput$ownTarget

# .. and also calculate the deviation to the consensus corridor:
dLowInput$corridorMedian <- apply(dLowInput,1,function(x) consensus$median[x["regionName"],x["sampleName"]] )
dLowInput$corridorLower <- apply(dLowInput,1,function(x) consensus$lower[x["regionName"],x["sampleName"]] )
dLowInput$corridorUpper <- apply(dLowInput,1,function(x) consensus$upper[x["regionName"],x["sampleName"]] )
dLowInput$corridorDeviation <- ifelse(dLowInput$value>=dLowInput$corridorLower & dLowInput$value<=dLowInput$corridorUpper,0,ifelse(dLowInput$value<dLowInput$corridorLower, dLowInput$value-dLowInput$corridorLower, dLowInput$value-dLowInput$corridorUpper))


### DEVIATION FROM CONSENSUS (Supplementary Figures 10.b and 11.b) ###

yRange <- ceiling(max(c(dLowInput$ownTargetDeviation,dLowInput$corridorDeviation),na.rm=TRUE))
for(curSeries in unique(dLowInput$series)) {
	p <- ggplot(dLowInput[dLowInput$series==curSeries,], aes(x=as.numeric(quantity), y=corridorDeviation, fill=factor(prettyLabel))) + geom_point() + stat_summary(fun.y = mean, geom = 'ribbon', fun.ymax = max, fun.ymin = min, .alpha = 0.05, alpha = 0.5) + defaultPlotTheme(flipX=TRUE) + xlab("DNA quantity [ng]") + ylab("Deviation from corridor") + scale_fill_manual(values=plotColLookup$prettyLabel, guide=FALSE) + scale_x_continuous(breaks=1:length(levels(dLowInput$quantity)),labels=levels(dLowInput$quantity)) + facet_grid(prettyLabel~., drop=FALSE) + geom_hline(y_intercept=0) + theme(axis.line.x=element_blank()) + ylim(-yRange,yRange)
	svgPlotGG(p, paste0("low_input_dev_consensus_",curSeries), 8,24)
}

### DEVIATION FROM OWN TARGET (Supplementary Figures 10.c and 11.c) ###

for(curSeries in unique(dLowInput$series)) {
	p <- ggplot(dLowInput[!is.na(dLowInput$ownTargetDeviation) & dLowInput$series==curSeries,], aes(x=as.numeric(quantity), y=ownTargetDeviation, fill=factor(prettyLabel))) + geom_point() + stat_summary(fun.y = mean, geom = 'ribbon', fun.ymax = max, fun.ymin = min, .alpha = 0.05, alpha = 0.5) + defaultPlotTheme(flipX=TRUE) + xlab("DNA quantity [ng]") + ylab("Deviation from own target value") + scale_fill_manual(values=plotColLookup$prettyLabel, guide=FALSE) + scale_x_continuous(breaks=1:length(levels(dLowInput$quantity)),labels=levels(dLowInput$quantity)) + facet_grid(prettyLabel~., drop=FALSE) + geom_hline(y_intercept=0) + theme(axis.line.x=element_blank()) + ylim(-yRange,yRange)
	svgPlotGG(p, paste0("low_input_dev_owntarget_",curSeries), 8,24)
}

### STATUS (Supplementary Figures 10.a and 11.a) ###

for(curSeries in unique(dLowInput$series)) {
	p <- ggplot(dLowInputSuccessRate[dLowInputSuccessRate$series==curSeries,], aes(x=quantity, y=successRate)) + geom_bar(stat="identity") + defaultPlotTheme(flipX=TRUE) + xlab("DNA quantity [ng]") + ylab("Success rate (successful / attempted measurements)")  + facet_grid(prettyLabel~., drop=FALSE) + geom_hline(y_intercept=0) + theme(axis.line.x=element_blank())
	svgPlotGG(p, paste0("low_input_dev_status_",curSeries), 5, 24)
}







###########################################################################
############ PROSTATE CANCER COHORT (Supplementary Figures 15) ############
###########################################################################

# prostate cancer cohort (and normal tissue samples) 
# 13_T01 ... 200_T01 	Prostate cancer samples
# XX-T01 ... XX-T06	Multiple sections of same tumor sample
# 001_N bis 008_N	Normal prostate tissue

# load in the assay data:
dCohort <- melt(lapply(lapply(c("EpiTyper"="EpiTyper_3","Infinium"="Infinium"), function(assayName) { read.table(paste0("data/results_prostate_cohort_", assayName, ".txt"), sep="\t", header=TRUE, check.names=FALSE, comment.char="", row.names=2) }), function(x) { colnames(x) <- toCamelCase(colnames(x)); x[x==""|x=="failed"]<-NA; x <- data.matrix(x[,-c(1,2,3)]); x[rowSums(!is.na(x))>0,]}))
colnames(dCohort) <- c("regionName","sampleName","methValue","assayName")

# every region and sample included in this table has been measured by both assays:
# sum(aggregate(assayName~regionName+sampleName,dCohort,length)[,3]<1, na.rm=T)

# extract sample and section numbers, as well as tumor-normal status from the sample names:
dCohort$sampleId <- gsub("^(\\d+[TN])(.+)$", "\\1", dCohort$sampleName, perl=TRUE)
dCohort$sectionId <- gsub("^(\\d+[TN])(.+)$", "\\2", dCohort$sampleName, perl=TRUE)
dCohort$`Sample type` <- ifelse(grepl("T",dCohort$sampleName), "Tumor", "Normal")

# assign colors:
tnColors <- structure(brewer.pal(3,"BrBG"),names=c("Tumor","Other","Normal"))
sectColors <- structure(brewer.pal(length(unique(dCohort$sectionId)),"Pastel1"),names=unique(dCohort$sectionId))

# plot scatter comparing the measurements of both assays,  ...
ggData <- dcast(dCohort, sampleName + regionName + `Sample type` + sampleId + sectionId ~ assayName, value.var="methValue")
ggData$regionName <- regionNameLookup[as.character(ggData$regionName)]
ggData2 <- sapply(split(ggData,ggData$regionName),function(x) { above <- sum(x$EpiTyper>x$Infinium,na.rm=T)>nrow(x)/2; c(
	"r"=paste("r =", round(cor(x$Infinium,x$EpiTyper,use="pairwise"),2)),
	"x"=quantile(range(x$Infinium,na.rm=TRUE),ifelse(!above,0.025,0.975)),
	"y"=quantile(range(x$EpiTyper,na.rm=TRUE),ifelse(!above,0.975,0.025)),
	"v"=ifelse(!above,1,0),
	"h"=ifelse(!above,0,1)
)})
ggData2 <- data.frame("regionName"=colnames(ggData2), t(ggData2))
ggData2$x <- as.numeric(ggData2$x)
ggData2$y <- as.numeric(ggData2$y)
ggData2$v <- as.numeric(ggData2$v)
ggData2$h <- as.numeric(ggData2$h)
# (a) ... per region
d <- ggplot(ggData, aes(x=Infinium*100, y=EpiTyper*100, color=`Sample type`)) + geom_point() + scale_color_manual(values=tnColors) + xlab(paste(DNA_METH_LABEL, "- Infinium")) + ylab(paste(DNA_METH_LABEL, "- EpiTyper 3")) + defaultPlotTheme(fontSize=9) + geom_text(aes(x=x*100,y=y*100,label=r,hjust=h,vjust=v),color="black",size=3,data=ggData2) + facet_wrap(~ regionName, scale="fixed",ncol=5)  + ylim(0,100) + xlim(0,100)
svgPlotGG(d, "cohort_scatter_per_region", 22, 9)
# (b) ... per sample
d <- ggplot(ggData, aes(x=Infinium*100, y=EpiTyper*100, color=sectionId)) + geom_point() + scale_color_manual(values=sectColors) + xlab(paste(DNA_METH_LABEL, "- Infinium")) + ylab(paste(DNA_METH_LABEL, "- EpiTyper 3")) + defaultPlotTheme() + facet_wrap(~ sampleId, scale="fixed",ncol=8)   + ylim(0,100) + xlim(0,100)
svgPlotGG(d, "cohort_scatter_per_sample", 26, 62)

# TRAIN CLASSIFIERS FOR DISCRIMINATING TUMOR/NORMAL SAMPLES:

nRndCombos <- 1000
cmp <- sapply(split(dCohort$sampleName, as.factor(dCohort$`Sample type`)), function(x) as.character(unique(x)))
nTrainSet <- c("Tumor"=8,"Normal"=4)
allSamples <- as.character(unlist(cmp))

# determine all possible combinations (n=70) of four out of eight normal samples:
allNormalCombos <- unlist(combn(cmp$Normal,nTrainSet["Normal"]))
# .. then pair each combination up with 1000 random selections of tumor samples to build training/test sets:
combos <- data.frame(matrix(NA, ncol=length(allSamples)-sum(nTrainSet), nrow=ncol(allNormalCombos)*nRndCombos), stringsAsFactors=FALSE)
for(i in 1:ncol(allNormalCombos)) {
	off <- ((i-1)*nRndCombos)
	combos[off+(1:nRndCombos),] <- t(sapply(1:nRndCombos, function(x) {
		setdiff(allSamples,c(allNormalCombos[,i],sample(cmp$Tumor,nTrainSet["Tumor"])))
	}))
}

# comment out for an alternative run with big training and small test sets:
#combos <- as.data.frame(t(apply(combos,1,function(x) setdiff(allCmp,x))))

# run the classifiers:
classifierResults <- list()
for(curAssay in unique(dCohort$assayName)) {
	message("\t*", curAssay)
	classifierData <- dcast(dCohort[dCohort$assayName==curAssay,], sampleName ~ regionName, value.var="methValue")
	classifierData$sampleName <- as.character(classifierData$sampleName)
	rownames(classifierData) <- classifierData$sampleName
	x <- data.frame("Target"=grepl("N",rownames(classifierData))*1,classifierData[,-1])

	# N.B. set nTrials=1, since we already determined the random trials by setting the combos and we're not adding any noise:
	classifierResults[[curAssay]] <- runClassifierTrials(x,nTrials=1,comparisons=cmp,testTrainCombos=combos)
}

# plot ROC curves (Supplementary Figures 15.b):
classifierResultsX <- classifierResults
perfPlotCols <- plotColLookup$assayGroup
svgPlot("cohort_roc", 3.5, 4)
for(n in names(classifierResultsX)) {
	plot(classifierResultsX[[n]]$perf, add=which(names(classifierResultsX)==n)!=1, col=perfPlotCols[n], bty="l", avg="threshold")
}
legend("bottomright", names(classifierResultsX), col=perfPlotCols[names(classifierResultsX)], lty=1, bty="n",lwd=2)
dev.off()

