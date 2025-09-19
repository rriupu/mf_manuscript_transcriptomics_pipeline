library(ComplexUpset)
library(tidyverse)
library(optparse)

option_list = list(

  make_option(c("-i", "--input_dir"), type = "character", default = NULL, 
              help = "Directory with the results of the TFBS enrichment analysis (Mandatory).", 
              metavar = "character"),

  make_option(c("-c", "--comparison"), type = "character", default = NULL, 
              help = "Comparison for which to generate the plots (Mandatory).", 
              metavar = "character"),
  
  make_option(c("-o", "--output_dir"), type = "character", default = NULL, 
              help = "Output directory (Mandatory). ", 
              metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

input_dir = opt$input_dir
comparison_to_extract = opt$comparison
output_dir = opt$output

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

# Prepare df for upset plot
# Rows: dataset, columns: TFs
all_tfs = unique(do.call(c, significant_tfs_ls))

upset_ls = lapply(significant_tfs_ls, function(x) {return(all_tfs %in% x)})
upset_df = as.data.frame(do.call(cbind, upset_ls))

comparisons = basename(dirname(dirname(files)))

colnames(upset_df) = comparisons
rownames(upset_df) = all_tfs

upset_df_vs_M1 = upset_df[, grepl("vs_M1", colnames(upset_df))]
to_rm = sapply(1:nrow(upset_df_vs_M1), function(x) {
    rm = all(upset_df_vs_M1[x, ] == FALSE)
    return(rm)
})
upset_df_vs_M1 = upset_df_vs_M1[!to_rm, ]

upset_df_vs_untreated = upset_df[, grepl("vs_Untreated", colnames(upset_df))]
to_rm = sapply(1:nrow(upset_df_vs_untreated), function(x) {
    rm = all(upset_df_vs_untreated[x, ] == FALSE)
    return(rm)
})
upset_df_vs_untreated = upset_df_vs_untreated[!to_rm, ]

# Plot
pdf(file.path(output_dir, "upset_plot_all_comparisons.pdf"), width = 18, height = 12)
upset(upset_df, intersect = colnames(upset_df), name = "Comparisons")
dev.off()

pdf(file.path(output_dir, "upset_plot_vs_M1.pdf"), width = 18, height = 12)
upset(upset_df_vs_M1, intersect = colnames(upset_df_vs_M1), name = "Comparisons")
dev.off()

pdf(file.path(output_dir, "upset_plot_vs_Untreated.pdf"), width = 18, height = 12)
upset(upset_df_vs_untreated, intersect = colnames(upset_df_vs_untreated), name = "Comparisons")
dev.off()