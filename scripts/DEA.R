## --------------- ##
## Library loading ##
## --------------- ##

library(optparse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(limma)
library(dplyr)
library(ggh4x)

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

dset$condition = relevel(dset$condition, ref = "Untreated")

## ----------- ##
## DEA testing ##
## ----------- ##

dds = DESeq(dset)

pdf(file.path(outputDir, "dispersion_estimates.pdf"))
plotDispEsts(dds)
dev.off()

conditions = unique(conditionMapping$condition)

# ---------------------- #
# -- Untreated vs all -- #
# ---------------------- #
conditions_contrast1 = conditions[!conditions %in% "Untreated"]
for (condition in conditions_contrast1) {

    message(paste0("Running comparison: ", condition, " vs Untreated"))
    run_comparison(
        deseq_dset = dds,
        conditions = c(condition, "Untreated"),
        log2FCthr = c(-0.58, 0.58),
        pvalThr = 0.05,
        outputDir = outputDir)

}

# --------------- #
# -- M1 vs all -- #
# --------------- #
dds$condition = relevel(dds$condition, ref = "M1")
dds = nbinomWaldTest(dds)
conditions_contrast2 = conditions[!conditions %in% "M1"]
for (condition in conditions_contrast2) {

    message(paste0("Running comparison: ", condition, " vs M1"))
    run_comparison(
        deseq_dset = dds,
        conditions = c(condition, "M1"),
        log2FCthr = c(-0.58, 0.58),
        pvalThr = 0.05,
        outputDir = outputDir)

}

# **** Create summary of volcano plots **** #

# Get all files
all_results = list.files(
    path = outputDir,
    pattern = "DEA_results.tsv",
    recursive = T,
    full.names = T
)

# Parse the files to create a table for plotting
dfs_volcano_ls = list()
dfs_heatmap_058_ls = list()
dfs_heatmap_15_ls = list()
dfs_heatmap_2_ls = list()

for (i in 1:length(all_results)) {

    path_parts = strsplit(all_results[i], "/", fixed = T)
    model = path_parts[[1]][3]
    comparison = path_parts[[1]][4]
    conditions = strsplit(comparison, "_vs_")
    condition1 = conditions[[1]][1]
    condition2 = conditions[[1]][2]


    df = read.delim(all_results[i], header = T) %>%
        select(gene, log2FoldChange, pvalue, padj) %>%
        mutate(model = model, comparison = comparison, condition1 = condition1, condition2 = condition2)
    # dfs_volcano_ls[[i]] = df 

    n_degs_058 = df %>%
        filter(abs(log2FoldChange) >= 0.58, pvalue <= 0.05) %>%
        select(model, condition1, condition2) %>%
        summarise(n_degs = n()) %>%
        mutate(model = model, condition1 = condition1, condition2 = condition2)

    dfs_heatmap_058_ls[[i]] = n_degs_058

    n_degs_15 = df %>%
        filter(abs(log2FoldChange) >= 1.5, pvalue <= 0.05) %>%
        select(model, condition1, condition2) %>%
        summarise(n_degs = n()) %>%
        mutate(model = model, condition1 = condition1, condition2 = condition2)

    dfs_heatmap_15_ls[[i]] = n_degs_15

    n_degs_2 = df %>%
        filter(abs(log2FoldChange) >= 2, pvalue <= 0.05) %>%
        select(model, condition1, condition2) %>%
        summarise(n_degs = n()) %>%
        mutate(model = model, condition1 = condition1, condition2 = condition2)

    dfs_heatmap_2_ls[[i]] = n_degs_2

}

# supervolcano_df = do.call(rbind, dfs_volcano_ls)
# supervolcano_plot = ggplot(supervolcano_df, aes(x = log2FoldChange, y = -log10(pvalue), size = )) + 
#     geom_point() +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
#     geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") + 
#     facet_nested(model + condition2 ~ condition1) + 
#     theme_classic()

# ggsave(file.path(outputDir, "supervolcano_plot.jpg"), supervolcano_plot)

heatmap_058_df = do.call(rbind, dfs_heatmap_058_ls)
degs_heatmap = ggplot(heatmap_058_df, aes(x = condition1, y = n_degs, label = n_degs)) + 
    geom_bar(stat = "identity") +
    geom_text(size = 2, vjust = -0.5, hjust = 0.5) +
    facet_grid(model ~ condition2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(file.path(outputDir, "contrasts_summary_barplot_lfc_058.pdf"), degs_heatmap, width = 10)

heatmap_15_df = do.call(rbind, dfs_heatmap_15_ls)
degs_heatmap = ggplot(heatmap_15_df, aes(x = condition1, y = n_degs, label = n_degs)) + 
    geom_bar(stat = "identity") +
    geom_text(size = 2, vjust = -0.5, hjust = 0.5) +
    facet_grid(model ~ condition2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(file.path(outputDir, "contrasts_summary_barplot_lfc_15.pdf"), degs_heatmap, width = 10)

heatmap_2_df = do.call(rbind, dfs_heatmap_2_ls)
degs_heatmap = ggplot(heatmap_2_df, aes(x = condition1, y = n_degs, label = n_degs)) + 
    geom_bar(stat = "identity") +
    geom_text(size = 2, vjust = -0.5, hjust = 0.5) +
    facet_grid(model ~ condition2) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(file.path(outputDir, "contrasts_summary_barplot_lfc_2.pdf"), degs_heatmap, width = 10)






dfs_heatmap_ls = list()

ls_index_counter = 1

for (i in 1:length(all_results)) {

    ls_index_counter = ls_index_counter + 1
    path_parts = strsplit(all_results[i], "/", fixed = T)
    model = path_parts[[1]][3]
    comparison = path_parts[[1]][4]
    conditions = strsplit(comparison, "_vs_")
    condition1 = conditions[[1]][1]
    condition2 = conditions[[1]][2]

    df = read.delim(all_results[i], header = T) %>%
        select(gene, log2FoldChange, pvalue, padj) %>%
        mutate(model = model, comparison = comparison, condition1 = condition1, condition2 = condition2)

    for (log2FCthr in c(0.58, 1.5, 2)) {

        for (pvalThr in c(0.05, 0.01, 0.001)) {

            new_df = df %>%
                filter(abs(log2FoldChange) >= log2FCthr, padj <= pvalThr) %>%
                select(model, condition1, condition2) %>%
                summarise(n_degs = n()) %>%
                mutate(model = model, condition1 = condition1, condition2 = condition2, log2FCthr = log2FCthr, pvalThr = pvalThr)

            dfs_heatmap_ls[[ls_index_counter]] = new_df
            ls_index_counter = ls_index_counter + 1

        }

    }

}

plot_df = do.call(rbind, dfs_heatmap_ls)

pdf(file.path(outputDir, "summary_barplots.pdf"), width = 10)
for (model in unique(plot_df$model)) {

    p = ggplot(plot_df[plot_df$model == model, ], aes(x = condition1, y = n_degs, label = n_degs)) + 
        geom_bar(stat = "identity") +
        geom_text(size = 2, vjust = -0.5, hjust = 0.5) +
        facet_nested(log2FCthr + pvalThr ~ condition2) +
        ggtitle(model) +
        ylim(c(0, max(plot_df$n_degs[plot_df$model == model]) + 500)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

    print(p)

}
dev.off()


ultimate_barplot = ggplot(plot_df, aes(x = condition1, y = n_degs, label = n_degs, fill = model)) + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(size = 2, vjust = -0.5, hjust = 0.5, position = position_dodge(width = 0.9)) +
        facet_nested(log2FCthr + pvalThr ~ condition2) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(file.path(outputDir, "ultimate_barplot.pdf"), ultimate_barplot, width = 30)


# Reactome enrichment over models, thresholds, and conditions

reactomeEnrichment = function(foreground_genes_df, universe_genes_df, condition) {

    reactome_enrichment = enrichPathway(
        gene = foreground_genes_df$entrez,
        universe = universe_genes_df$entrez,
        organism = "human",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05,
        readable = T)

    if (nrow(fortify(reactome_enrichment, showCategory = 10, split = NULL)) > 0 && !is.null(reactome_enrichment)) {

        enrichment_result = reactome_enrichment %>%
            as.data.frame() %>%
            filter(p.adjust <= 0.05) %>%
            select(Description, pvalue, p.adjust) %>%
            rename(c("enrichment_pval" = "pvalue", "enrichment_padj" = "p.adjust"))

    } else {

        enrichment_result = data.frame(Description = character(), enrichment_pval = numeric(), enrichment_padj = numeric())

    }

    return(enrichment_result)

}

KEGGenrichment = function(foreground_genes_df, universe_genes_df, condition) {

    kegg_enrichment = enrichKEGG(
        gene = foreground_genes_df$entrez,
        universe = universe_genes_df$entrez,
        organism = 'hsa',
        pAdjustMethod = "BH",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05)

    if (nrow(fortify(kegg_enrichment, showCategory = 10, split = NULL)) > 0 && !is.null(kegg_enrichment)) {
        
        kegg_res_df = as.data.frame(kegg_enrichment) %>%
            select(Description, pvalue, p.adjust) %>%
            rename(c("enrichment_pval" = "pvalue", "enrichment_padj" = "p.adjust"))

    } else {

        kegg_res_df = data.frame(Description = character(), enrichment_pval = numeric(), enrichment_padj = numeric())

    }

    return(kegg_res_df)

}


reactome_enrich_ls = list()
kegg_enrich_ls = list()
ls_index_counter = 1

for (i in 1:length(all_results)) {

    print(paste0(i, "/", length(all_results)))
    path_parts = strsplit(all_results[i], "/", fixed = T)
    model = path_parts[[1]][3]
    comparison = path_parts[[1]][4]
    conditions = strsplit(comparison, "_vs_")
    condition1 = conditions[[1]][1]
    condition2 = conditions[[1]][2]

    df = read.delim(all_results[i], header = T) %>%
        select(gene, log2FoldChange, pvalue, padj) %>%
        mutate(model = model, comparison = comparison, condition1 = condition1, condition2 = condition2)

    universe_genes = mapGeneSymbolsToEntrez(df$gene)

    for (log2FCthr in c(0.58, 1.5, 2)) {

        for (pvalThr in c(0.05, 0.01, 0.001)) {

            new_df = df %>%
                filter(log2FoldChange >= log2FCthr, padj <= pvalThr)

            condition_genes = new_df %>% pull(gene)
            if (length(condition_genes) > 0) {

                condition_genes = mapGeneSymbolsToEntrez(condition_genes)
                enrichment_results = reactomeEnrichment(condition_genes, universe_genes, condition1)
                kegg_enrichment_results = KEGGenrichment(condition_genes, universe_genes, condition1)

            } else {

                enrichment_results = data.frame(Description = character(), enrichment_pval = numeric(), enrichment_padj = numeric())
                kegg_enrichment_results = data.frame(Description = character(), enrichment_pval = numeric(), enrichment_padj = numeric())

            }

            enrichment_results = enrichment_results %>%
                mutate(model = model, condition1 = condition1, condition2 = condition2, log2FCthr = log2FCthr, pvalThr = pvalThr)

            kegg_enrichment_results = kegg_enrichment_results %>%
                mutate(model = model, condition1 = condition1, condition2 = condition2, log2FCthr = log2FCthr, pvalThr = pvalThr)

            reactome_enrich_ls[[ls_index_counter]] = enrichment_results
            kegg_enrich_ls[[ls_index_counter]] = kegg_enrichment_results
            ls_index_counter = ls_index_counter + 1

        }

    }

}

plot_df = do.call(rbind, reactome_enrich_ls)

pathways = unique(plot_df$Description)
plot_df$y = 1
pdf(file.path(outputDir, "reactome_enrichments_4_donors_with_donor_effect_vs_Untreated_log2FC_1_pval_005.pdf"))
for (pathway in pathways) {

    p = ggplot(plot_df[plot_df$Description == pathway, ], aes(x = condition1, y = y)) +
    geom_bar(stat = "identity", position = "dodge") +
    # facet_nested(log2FCthr + pvalThr ~ condition2) +
    ggtitle(pathway) +
    ylab("") +
    xlab("") +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 6
        ),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
    ) + 
    coord_equal()

    print(p)

}
dev.off()


# goldilocks plot

final_df = data.frame(
    Description = character(),
    enrichment_padj = numeric(),
    condition1 = character(),
    y = numeric()
)

descriptions = unique(plot_df$Description)
conditions = unique(conditionMapping$condition)
conditions = conditions[!conditions %in% "Untreated"]

for (description in descriptions) {

    rows_with_data = plot_df %>% filter(Description == description) %>% select(Description, enrichment_padj, condition1, y)
    conditions_with_data = rows_with_data %>% pull(condition1)
    conditions_without_data = conditions[!(conditions %in% conditions_with_data)]

    rows_without_data = data.frame(
        Description = description,
        enrichment_padj = NA,
        condition1 = conditions_without_data,
        y = 0
    )

    final_df = rbind(final_df, rows_with_data, rows_without_data)

}

pathways = unique(final_df$Description)
pdf(file.path(outputDir, "kegg_goldilocks.pdf"))
for (pathway in pathways) {

    p = ggplot(final_df[final_df$Description == pathway, ], aes(x = condition1, y = y)) +
    geom_bar(stat = "identity", position = "dodge") +
    # facet_nested(log2FCthr + pvalThr ~ condition2) +
    ggtitle(pathway) +
    ylab("") +
    xlab("") +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 6
        ),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
    ) + 
    coord_equal()

    print(p)

}
dev.off()









ggplot(plot_df[plot_df$Description == pathways[62], ], aes(x = condition1, y = y, fill = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_nested(log2FCthr + pvalThr ~ condition2) +
    ggtitle(pathways[62]) +
    ylab("") +
    xlab("") +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 6
        ),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
    ) 
