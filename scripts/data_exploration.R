## --------------- ##
## Library loading ##
## --------------- ##

library(optparse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(limma)

source("scripts/utils.R")

## -------------- ##
## Read arguments ##
## -------------- ##

option_list = list(

  make_option(c("-c", "--conditionMapping"), type = "character", default = NULL, 
              help = "Text file containing a mapping between a sample ID and its group (Mandatory).", 
              metavar = "character"),
  
  make_option(c("-g", "--geneCountsDir"), type = "character", default = NULL, 
              help = "Directory containing the results of htseq-count (Mandatory). ", 
              metavar = "character"),
  
  make_option(c("-o", "--outputDir"), type = "character", default = NULL, 
              help = "Directory where results will be saved (Mandatory). ", 
              metavar = "character")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

## ---------------- ##
## Read input files ##
## ---------------- ##

conditionMapping = opt$conditionMapping
conditionMapping = prepareConditionMapping(conditionMapping)

geneCountsDir = opt$geneCountsDir
outputDir = opt$outputDir

## -------------------- ##
## Prepare DESeqDataset ##
## -------------------- ##

dset = prepareDESeqDataset(conditionMapping, geneCountsDir, formula = "~ donor + condition")

dset$condition = relevel(dset$condition, ref = "Untreated")

## ---------------- ##
## Data exploration ##
## ---------------- ##

# Variance stabilization
vsd = vst(dset, blind = TRUE)

# Batch effect removal
assay(vsd) = limma::removeBatchEffect(assay(vsd), vsd$donor)

# Sample-sample distance
plotSampleDistHeatmap(vsd, outputPath = file.path(outputDir, "sample_sample_distance_blind_no_donor_effect.pdf"))

# Sample-sample correlation
vsd_cor = cor(assay(vsd))
heatmap_cor = pheatmap(vsd_cor)
ggsave(file.path(outputDir, "correlation_heatmap_no_donor_effect.pdf"), heatmap_cor)

# PCA
pca_plot = plotPCA(vsd, c("condition", "donor"))
p = ggplot(pca_plot$data, aes(x = PC1, y = PC2, color = condition, shape = donor)) + geom_point()
ggsave(file.path(outputDir, "PCA_blind_no_donor_effect.pdf"), p, width = 14)

p2 = ggplot(pca_plot$data, aes(x = PC1, y = PC2, color = condition, shape = donor)) + geom_point() + facet_wrap(. ~ condition)
ggsave(file.path(outputDir, "PCA_blind_no_donor_effect_facet.pdf"), p2, width = 14)