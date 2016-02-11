#### Load Utility Functions ####
source("analysis/methBench_lib.R")

#### Project-Specific Constants ####
PROJECT <- "methBench"
ROOT <- getwd()
DEFAULT_LIBS <- c("RColorBrewer","gplots","fastcluster","ggplot2","reshape2","ROCR","e1071","randomForest","pheatmap","beeswarm","vcd","DMwR","Cairo")
for(lib in DEFAULT_LIBS) library(lib,character.only=TRUE)
REMOVE_REGIONS <- c("region_11") # region_11 contains a SNP
HIGH_COR_THRESHOLD <- 0.8
DNA_METH_LABEL <- "DNA methylation (%)"
MISSING_VALUES <- c("#DIV/0!"," ","","nd","NA","N/A","#N/A","ND","n.a.","Homozygous CpA","no PCR product")
CORRIDOR_EXTENSION <- 10
CORRIDOR_MIN_NUM_ASSAY_TYPES <- 3
CLASSIFIER_TRIALS <- 1000
NOISE_TYPES <- list("Random error"=addRandomError, "Uniform noise"=addUniformNoise) 
NOISE_LEVELS <- c(0.0,0.25,0.5,1.0)
NOISE_COLORS <- list("Random error"=c("black",rep(brewer.pal(3,"PuOr")[1],length(NOISE_LEVELS)-1)),"Uniform noise"=c("black",rep(brewer.pal(3,"PuOr")[3],length(NOISE_LEVELS)-1)))
STANDARD_NOISE <- 0.25
STANDARD_NOISE_TYPE <- "Uniform noise"
REGION_SELECTION_METHODS <- list("Random"=pickRandomRegions, "Top"=pickTopRegions)
REGION_NUMBERS <- c(5,3,1)
latestPrepName <- paste("methBench","afterPrep","latest","session.Rdata",sep="_")

#### Initialise Workspace ####
setwd(ROOT)
options(stringsAsFactors=FALSE)
options(max.print=250)
set.seed(1234) # set random seed to ensure reproducibility of results

#### Load Data ####

source("analysis/methBench_data.R")



# set a save point with all project-relevant data loaded:
save.image(paste("methBench","afterPrep",format(Sys.time(), "%Y%m%d_%H%M"),"session.Rdata",sep="_"))
save.image(latestPrepName)
save.image()

