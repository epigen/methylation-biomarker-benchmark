# Entry point script for teh analysis presented in the manuscript:
# "Quantitative comparison of locus-specific DNA methylation assays for large-scale validation studies and epigenetic biomarker development"

# We executed the code with R version 3.1.2 and the following 
# library / package versions relevant to the results:
# beeswarm_0.2.0
# Cairo_1.5-6
# DMwR_0.4.1
# e1071_1.6-4
# fastcluster_1.1.16
# ggplot2_1.0.1
# gplots_2.17.0
# lattice_0.20-30
# pheatmap_1.0.2
# randomForest_4.6-10
# reshape2_1.4.1
# RColorBrewer_1.1-2
# ROCR_1.0-7
# vcd_1.3-2

# load data, set constants, etc.:
source("analysis/methBench_prepare.R")

# generate plots and summary statistics linked to Figure 1-5:
source("analysis/methBench_fig1.R")
source("analysis/methBench_fig2.R")
source("analysis/methBench_fig3.R")
source("analysis/methBench_fig4.R")
source("analysis/methBench_fig5.R")

# generate plots and summary statistics for the supplementary section:
source("analysis/methBench_suppFigs.R")

# run additional analysis supplementary to the primary benchmark:
source("analysis/methBench_additions.R")



