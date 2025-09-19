library(optparse)
library(tidyverse)

option_list = list(

  make_option(c("-i", "--input_dir"), type = "character", default = NULL, 
              help = "Directory with the results of the TFBS enrichment analysis (Mandatory).", 
              metavar = "character"),

  make_option(c("-c", "--comparison"), type = "character", default = NULL, 
              help = "Comparison for which to generate the plots (Mandatory).", 
              metavar = "character"),
  
  make_option(c("-o", "--output_dir"), type = "character", default = NULL, 
              help = "Output plot (Mandatory). ", 
              metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

input_dir = opt$input_dir
comparison_to_extract = opt$comparison
output_dir = opt$output_dir

dir.create(output_dir, recursive = T, showWarnings = F)

get_significant_tfs = function(df, log10_pval_thr = 1.33) {

    sig_tfs = df %>%
        filter(pValueLog >= log10_pval_thr) %>%
        pull(collection) %>%
        unique()
    
    return(sig_tfs)

}

# Find all files
files = list.files(input_dir, pattern = "allEnrichments.tsv", full.names = T, recursive = T)
if (comparison_to_extract == "all") {

    files = files[!grepl("upregulated", files)]
    files = files[!grepl("downregulated", files)]

} else if (comparison_to_extract == "upregulated") {

    files = files[!grepl("/all/", files)]
    files = files[!grepl("downregulated", files)]

} else if (comparison_to_extract == "downregulated") {

    files = files[!grepl("/all/", files)]
    files = files[!grepl("upregulated", files)]

}

# Read into list
dfs_ls = lapply(files, read.delim)

# Get significant TFs
significant_tfs_ls = lapply(dfs_ls, get_significant_tfs)

all_tfs = unique(do.call(c, significant_tfs_ls))
comparisons = basename(dirname(dirname(files)))

# Prepare df for upset plot
# Rows: dataset, columns: TFs

plotting_ls = lapply(1:length(significant_tfs_ls), function(x) {
    df = data.frame(
        comparison = comparisons[x],
        tf = all_tfs,
        value = all_tfs %in% significant_tfs_ls[[x]]
    )
    return(df)
})

# Plot for all

plotting_df = as.data.frame(do.call(rbind, plotting_ls))
plotting_df$tf = factor(plotting_df$tf, levels = sort(all_tfs, decreasing = F))
plotting_df$comparison = factor(plotting_df$comparison, levels = sort(comparisons, decreasing = T))

plotting_df$max_pval = sapply(1:nrow(plotting_df), function(x) {

    comparison = plotting_df$comparison[x]
    tf = plotting_df$tf[x]

    max_pval = dfs_ls[[which(comparisons %in% comparison)]] %>%
        filter(collection == tf) %>%
        pull(pValueLog) %>%
        max()

    return(max_pval)

})

# Get max and min pvals from the significant ones
signif_pvals = plotting_df$max_pval[plotting_df$value == TRUE]
max_pval = max(signif_pvals)
min_pval = min(signif_pvals)

# Conditional minmax transformation
plotting_df$max_pval_scaled = NA
plotting_df$max_pval_scaled[plotting_df$value == TRUE] = (plotting_df$max_pval[plotting_df$value == TRUE] - min_pval) / (max_pval - min_pval)
plotting_df$max_pval_scaled[plotting_df$value == FALSE] = 0
plotting_df$alpha = ifelse(plotting_df$value == TRUE, yes = 1, no = 0)

dotplot_all = ggplot(plotting_df, aes(x = tf, y = comparison, size = max_pval_scaled, alpha = alpha, color = value)) + 
    geom_point() +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 6
        ),
        legend.position = "none"
    )

ggsave(file.path(output_dir, "summary_dotplot_all.pdf"), dotplot_all, width = 18, height = 10)

p_all = ggplot(plotting_df, aes(x = tf, y = comparison, fill = value)) + 
    geom_tile(aes(width = 0.8, height = 0.8)) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    coord_fixed() +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 6
        ),
        legend.position = "none"
    )

ggsave(file.path(output_dir, "summary_plot_all.pdf"), p_all, width = 14)

# Plot for comparisons against M1

plotting_df_vs_M1 = plotting_df %>%
    filter(grepl("_vs_M1", comparison))
tfs_to_rm = sapply(all_tfs, function(x) {
    rm = all(plotting_df_vs_M1$value[plotting_df_vs_M1$tf == x] == FALSE)
    return(rm)
})
plotting_df_vs_M1 = plotting_df_vs_M1[!tfs_to_rm, ]

dotplot_M1 = ggplot(plotting_df_vs_M1, aes(x = tf, y = comparison, size = max_pval_scaled, alpha = alpha, color = value)) + 
    geom_point() +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 6
        ),
        legend.position = "none"
    )

ggsave(file.path(output_dir, "summary_dotplot_vs_M1.pdf"), dotplot_M1, width = 18, height = 10)

p_M1 = ggplot(plotting_df_vs_M1, aes(x = tf, y = comparison, fill = value)) + 
    geom_tile(aes(width = 0.8, height = 0.8)) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    coord_fixed() +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 6
        ),
        legend.position = "none"
    )

ggsave(file.path(output_dir, "summary_plot_vs_M1.pdf"), p_M1, width = 14)

# Plot for comparisons against Untreated

plotting_df_vs_Untreated = plotting_df %>%
    filter(grepl("_vs_Untreated", comparison))
tfs_to_rm = sapply(all_tfs, function(x) {
    rm = all(plotting_df_vs_Untreated$value[plotting_df_vs_Untreated$tf == x] == FALSE)
    return(rm)
})
plotting_df_vs_Untreated = plotting_df_vs_Untreated[!tfs_to_rm, ]

dotplot_untreated = ggplot(plotting_df_vs_Untreated, aes(x = tf, y = comparison, size = max_pval_scaled, alpha = alpha, color = value)) + 
    geom_point() +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 6
        ),
        legend.position = "none"
    )

ggsave(file.path(output_dir, "summary_dotplot_vs_Untreated.pdf"), dotplot_untreated, width = 18, height = 10)

p_untreated = ggplot(plotting_df_vs_Untreated, aes(x = tf, y = comparison, fill = value)) + 
    geom_tile(aes(width = 0.8, height = 0.8)) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
    coord_fixed() +
    xlab("") +
    ylab("") +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 6
        ),
        legend.position = "none"
    )

ggsave(file.path(output_dir, "summary_plot_vs_Untreated.pdf"), p_untreated, width = 14)