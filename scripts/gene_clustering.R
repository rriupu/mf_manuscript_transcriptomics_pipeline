## --------------- ##
## Library loading ##
## --------------- ##

library(data.table)
library(tidyverse)
library(pheatmap)
library(optparse)
library(ComplexHeatmap)

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

dir.create(outputDir, recursive = T, showWarnings = F)

## -------------------- ##
## Prepare DESeqDataset ##
## -------------------- ##

dset = prepareDESeqDataset(conditionMapping, geneCountsDir, formula = "~ donor + condition")

smallestGroupSize = 4
keep = rowSums(counts(dset) >= 10) >= smallestGroupSize
dset = dset[keep, ]

## -------------------------- ##
## Subset columns of interest ##
## -------------------------- ##

cols_of_int = c(
    grep("_TCM", colnames(dset)),
    grep("_K_", colnames(dset)),
    grep("_TGFb", colnames(dset)),
    grep("_M1_", colnames(dset)),
    grep("_M2", colnames(dset)),
    grep("_ADO_", colnames(dset))
)

dset = dset[, unique(cols_of_int)]

## -------------------- ##
## Batch effect removal ##
## -------------------- ##

n_top_genes = 1000

vsd = vst(dset, blind = TRUE)
assay(vsd) = limma::removeBatchEffect(assay(vsd), vsd$donor)

## --------------------------------------------------------------------- ##
## Use top 1000 most variable genes across conditions for the clustering ##
## --------------------------------------------------------------------- ##

row_var = rowVars(assay(vsd))

top_genes = names(sort(row_var, decreasing = T)[1:n_top_genes])
vsd_top_genes = vsd[top_genes, ]

vsd_top_genes = assay(vsd_top_genes)

sample_annots = data.frame(
    row.names = colnames(vsd_top_genes),
    sample = sapply(strsplit(colnames(vsd_top_genes), "_"), "[[", 2),
    donor = sapply(strsplit(colnames(vsd_top_genes), "_"), "[[", 1)
)

# Z-scoring
row_means = rowMeans(vsd_top_genes)
row_vars = rowVars(vsd_top_genes)
vsd_top_genes_zscore = (vsd_top_genes - row_means)/row_vars

## ----------------- ##
## Building heatmaps ##
## ----------------- ##

gene_clusters = 4
column_clusters = 9

heatmap = Heatmap(
    vsd_top_genes_zscore,
    row_km = gene_clusters,
    column_km = column_clusters,
    top_annotation = HeatmapAnnotation(
        Condition = sample_annots$sample,
        Donor = sample_annots$donor
    ),
    show_row_names = F,
    column_names_gp = grid::gpar(fontsize = 8)
)

heatmap = draw(heatmap)

rd = row_order(heatmap)
cd = column_order(heatmap)

## ---------------------------------------- ##
## Run cluster-specific pathway enrichments ##
## ---------------------------------------- ##

cluster_enrich_ls = list()
for (cluster in 1:length(rd)) {

    rd_mask = as.numeric(names(rd)) == cluster
    rds = unlist(rd[rd_mask])

    foreground_genes_df =  mapGeneSymbolsToEntrez(rownames(vsd_top_genes_zscore)[rds])
    universe_genes_df = mapGeneSymbolsToEntrez(rownames(assay(vsd)))

    enrichment_results = enrichPathway(
        gene = foreground_genes_df$entrez,
        universe = universe_genes_df$entrez,
        organism = "human",
        readable = T)@result %>%
        filter(p.adjust <= 0.1) %>%
        select(Description, pvalue, p.adjust, qvalue, GeneRatio, BgRatio)
        
    cluster_enrich_ls[[cluster]] = enrichment_results

    out_dir = file.path(outputDir, cluster)
    dir.create(out_dir, recursive = T)
    
    genes = data.frame(
        gene = rownames(vsd_top_genes_zscore)[rds]
    )

    write.table(
        genes,
        file.path(out_dir, "genes.txt"),
        sep = "\t",
        col.names = T,
        row.names = F,
        quote = F
    )

    write.table(
        enrichment_results,
        file.path(out_dir, "enriched_pathways.txt"),
        sep = "\t",
        col.names = T,
        row.names = F,
        quote = F
    )

}

saveRDS(cluster_enrich_ls, file.path(outputDir, "cluster_enrichment_ls.rds"))

pdf(file.path(outputDir, "clusters.pdf"), width = 16)
print(heatmap)
dev.off()
