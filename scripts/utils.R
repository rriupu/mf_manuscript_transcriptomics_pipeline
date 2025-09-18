library(org.Hs.eg.db)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)
library(pheatmap)
library(RColorBrewer)

prepareConditionMapping = function(path) {

    conditionMapping = read.table(path, header = T)
    conditionMapping$donor = sapply(strsplit(conditionMapping$sample_id, "_"), "[[", 1)
    conditionMapping$condition = gsub("-", "_", conditionMapping$condition)

    return(conditionMapping)

}

prepareDESeqDataset = function(conditionMapping, geneCountsDir, formula) {

    geneCountFiles = list.files(geneCountsDir, pattern = "_geneCounts.tsv")

    sample_idx = sapply(conditionMapping$sample_id, function(x) {
        grep(x, geneCountFiles)
    })

    sampleTable = data.frame(
        sampleName = conditionMapping$sample_id,
        fileName = geneCountFiles[sample_idx],
        condition = factor(conditionMapping$condition),
        donor = factor(conditionMapping$donor))

    dset = DESeqDataSetFromHTSeqCount(
        sampleTable = sampleTable,
        directory = geneCountsDir,
        # design = ~ donor + condition)
        design = formula(formula))

    return(dset)

}

plotNReadsPerSample = function(dseq_dset, outputDir) {

    n_reads_per_sample = colSums(assay(dseq_dset))
    plotting_df = data.frame(sample = names(n_reads_per_sample), n_reads = n_reads_per_sample)
    n_reads_per_sample_barplot = ggplot(plotting_df, aes(x = sample, y = n_reads)) + 
        geom_bar(stat = "identity") + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    n_reads_per_sample_violin = ggplot(plotting_df, aes(x = 0, y = n_reads)) + 
        geom_violin() + 
        geom_dotplot(binaxis= "y", stackdir = "center", dotsize = 0.5, binwidth = (max(n_reads_per_sample) - min(n_reads_per_sample)) / 50) +
        xlab("") +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

    ggsave(file.path(outputDir, "n_reads_per_sample_barplot.pdf"), n_reads_per_sample_barplot, width = 14)
    ggsave(file.path(outputDir, "n_reads_per_sample_violin.pdf"), n_reads_per_sample_violin)

}

plotSampleDistHeatmap = function(vsd, outputPath) {

    sampleDists = dist(t(assay(vsd)))
    sampleDistMatrix = as.matrix(sampleDists)
    rownames(sampleDistMatrix) = paste0(vsd$donor, "_", vsd$condition)
    colnames(sampleDistMatrix) = NULL
    colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
    htmp = pheatmap(sampleDistMatrix,
            clustering_distance_rows = sampleDists,
            clustering_distance_cols = sampleDists,
            col = colors)

    ggsave(outputPath, htmp, height = 14)

}

plotVolcano = function(df, log2FCthr = c(-2, 2), pvalThr = 0.05, outputPath) {

    p = ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(size = 0.3) +
        geom_vline(xintercept = log2FCthr, linetype = "dashed") +
        geom_hline(yintercept = -log10(pvalThr), linetype = "dashed") +
        theme_classic()

    ggsave(outputPath, p)

}

mapGeneSymbolsToEntrez = function(genes) {

    map = data.frame(symbol = genes) %>%
        mutate(
            entrez = mapIds(org.Hs.eg.db, keys = genes,
            column = "ENTREZID",
            keytype = "SYMBOL")) %>% 
        filter(!is.na(entrez))
    
    return(map)

}

runGOenrichment = function(foreground_genes_df, universe_genes_df, condition, outputDir) {

    go_enrichment = enrichGO(
            gene = foreground_genes_df$entrez,
            universe = universe_genes_df$entrez,
            OrgDb = org.Hs.eg.db,
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.01,
            qvalueCutoff = 0.05,
            readable = TRUE)

    if (nrow(fortify(go_enrichment, showCategory = 10, split = NULL)) > 0 && !is.null(go_enrichment)) {
        dotplot_go_enrichment = dotplot(go_enrichment)
        ggsave(file.path(outputDir, paste0(condition, "_GO_enrichment.pdf")), dotplot_go_enrichment)

        enrichment_result = go_enrichment %>%
            filter(p.adjust <= 0.05)

        write.table(
            enrichment_result,
            file.path(outputDir, paste0(condition, "_GO_enrichment.tsv")),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")

    }

}

runGSEA = function(results_df, universe_genes_df, condition, outputDir) {

    annotated_genes = merge(results_df, universe_genes_df, by.x = "gene", by.y = "symbol")
    gene_list = annotated_genes$log2FoldChange
    names(gene_list) = annotated_genes$entrez
    gene_list = sort(gene_list, decreasing = T)

    gsea = gseGO(
        geneList = gene_list,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pvalueCutoff = 0.01)

    if (nrow(gsea@result) > 0 && !is.null(gsea)) {
        dotplot_gsea = dotplot(gsea)
        ggsave(file.path(outputDir, paste0(condition, "_GSEA.pdf")), dotplot_gsea)

        enrichment_result = gsea %>%
            filter(p.adjust <= 0.05)

        write.table(
            enrichment_result,
            file.path(outputDir, paste0(condition, "_GSEA.tsv")),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")

    }

}

runKEGGenrichment = function(foreground_genes_df, universe_genes_df, condition, outputDir) {

    kegg_enrichment = enrichKEGG(
        gene = foreground_genes_df$entrez,
        universe = universe_genes_df$entrez,
        organism = 'hsa',
        pAdjustMethod = "BH",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05)

    if (nrow(fortify(kegg_enrichment, showCategory = 10, split = NULL)) > 0 && !is.null(kegg_enrichment)) {
        kegg_res_df = as.data.frame(kegg_enrichment) %>%
            select(Description, GeneRatio, p.adjust, Count) %>%
            filter(p.adjust < 0.05)

        kegg_res_df$GeneRatio = sapply(kegg_res_df$GeneRatio, function(x) {eval(parse(text=x))})
        kegg_res_df = kegg_res_df[order(kegg_res_df$GeneRatio, decreasing = F), ]
        kegg_res_df$Description = factor(kegg_res_df$Description, levels = kegg_res_df$Description)

        kegg_dotplot = ggplot(kegg_res_df, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
            geom_point() +
            scale_color_continuous() +
            theme_light()

        ggsave(file.path(outputDir, paste0(condition, "_kegg_enrichment.pdf")), kegg_dotplot)

        enrichment_result = kegg_enrichment %>%
            filter(p.adjust <= 0.05)

        write.table(
            enrichment_result,
            file.path(outputDir, paste0(condition, "_KEGG_enrichment.tsv")),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")

    }

}

runReactomeEnrichment = function(foreground_genes_df, universe_genes_df, condition, outputDir) {

    reactome_enrichment = enrichPathway(
        gene = foreground_genes_df$entrez,
        universe = universe_genes_df$entrez,
        organism = "human",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05,
        readable = T)

    if (nrow(fortify(reactome_enrichment, showCategory = 10, split = NULL)) > 0 && !is.null(reactome_enrichment)) {

        dotplot_reactome_enrichment = dotplot(reactome_enrichment)
        ggsave(file.path(outputDir, paste0(condition, "_reactome_enrichment.pdf")), dotplot_reactome_enrichment)

        enrichment_result = reactome_enrichment %>%
            filter(p.adjust <= 0.05)

        write.table(
            enrichment_result,
            file.path(outputDir, paste0(condition, "_Reactome_enrichment.tsv")),
            col.names = T,
            row.names = T,
            quote = F,
            sep = "\t")

    }

}

runEnrichmentAnalyses = function(df, log2FCthr, pvalThr, conditions, outputDir) {

    condition1_genes = df$gene[df$padj <= pvalThr & df$log2FoldChange >= log2FCthr[2]]
    condition2_genes = df$gene[df$padj <= pvalThr & df$log2FoldChange <= log2FCthr[1]]

    universe_genes = mapGeneSymbolsToEntrez(df$gene)

    if (length(condition1_genes) > 0) {

        condition1_genes = mapGeneSymbolsToEntrez(condition1_genes)
        runGOenrichment(condition1_genes, universe_genes, conditions[1], outputDir)
        runKEGGenrichment(condition1_genes, universe_genes, conditions[1], outputDir)
        runReactomeEnrichment(condition1_genes, universe_genes, conditions[1], outputDir)
        runGSEA(df, universe_genes, conditions[1], outputDir)

        # # TODO add TFBS

    }

    if (length(condition2_genes) > 0) {

        condition2_genes = mapGeneSymbolsToEntrez(condition2_genes)
        runGOenrichment(condition2_genes, universe_genes, conditions[2], outputDir)
        runKEGGenrichment(condition2_genes, universe_genes, conditions[2], outputDir)
        runReactomeEnrichment(condition2_genes, universe_genes, conditions[2], outputDir)
        
    }

}

run_comparison = function(deseq_dset, conditions, log2FCthr = c(-2, 2), pvalThr = 0.05, outputDir) {

    outputDir = file.path(
        outputDir,
        paste0(conditions[1], "_vs_",
        conditions[2]))

    dir.create(
        outputDir,
        recursive = T,
        showWarnings = F)

    res_unshrunken = results(deseq_dset, contrast = c("condition", conditions[1], conditions[2]), alpha = 0.05)
    res_shrunken = lfcShrink(deseq_dset, coef = paste0("condition_", conditions[1], "_vs_", conditions[2]))

    pdf(file.path(outputDir, "MA_plots.pdf"))
    DESeq2::plotMA(res_unshrunken)
    DESeq2::plotMA(res_shrunken)
    dev.off()

    df = na.omit(res_shrunken) %>%
        as.data.frame() %>%
        rownames_to_column("gene")

    # sigOE = df %>%
    #     data.frame() %>%
    #     rownames_to_column(var = "gene") %>%
    #     filter(padj < pvalThr & abs(log2FoldChange) >= max(log2FCthr))

    # normalized_counts = counts(deseq_dset, normalized=T) %>% 
    #     data.frame() %>%
    #     rownames_to_column(var="gene") %>%
    #     as_tibble()
    
    # norm_OEsig = normalized_counts %>% 
    #           filter(gene %in% sigOE$gene)  

    # heat_colors <- brewer.pal(6, "YlOrRd")
    # sig_gene_heatmap = pheatmap(norm_OEsig[2:ncol(norm_OEsig)], 
    #     color = heat_colors, 
    #     cluster_rows = T, 
    #     show_rownames = F,
    #     border_color = NA, 
    #     fontsize = 10, 
    #     scale = "row", 
    #     fontsize_row = 10, 
    #     height = 20)
    # ggsave(file.path(outputDir, "significant_genes_heatmap.pdf"), sig_gene_heatmap)
    

    plotVolcano(df, log2FCthr = log2FCthr, pvalThr = pvalThr, outputPath = file.path(outputDir, "volcano_plot.pdf"))

    write.table(
        df,
        file.path(outputDir, "DEA_results.tsv"),
        col.names = T,
        row.names = F,
        quote = F,
        sep = "\t")

    runEnrichmentAnalyses(df, log2FCthr, pvalThr, conditions, outputDir)

}

plot_gene_boxplot = function(deseq_dset, gene_of_interest, conditions) {

    tcounts = as.data.frame(
        t((
            log2(
                counts(deseq_dset[gene_of_interest, ], normalized = TRUE, replaced = FALSE) +.5
            )
        ))
    ) %>%
    rownames_to_column("ID")

    colnames(tcounts)[2] = "expression"
    
    tcounts$condition = sapply(
        strsplit(tcounts$ID, "_"),
        "[[",
        2
    )

    tcounts = tcounts %>%
        filter(condition %in% c(conditions))

    ggplot(tcounts, aes(x = condition, y = expression)) + geom_boxplot() + ggtitle(gene_of_interest)

}