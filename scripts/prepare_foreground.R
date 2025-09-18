library(dplyr)
library(optparse)

## -------------- ##
## Read arguments ##
## -------------- ##

option_list = list(

  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "TSV file with one of the DEA results (Mandatory).", 
              metavar = "character"),

  make_option(c("-b", "--background"), type = "character", default = NULL, 
              help = "BED file with all the background promoter coordinates (Mandatory).", 
              metavar = "character"),

  make_option(c("-m", "--mapping"), type = "character", default = NULL, 
              help = "TSV file with an entrez ID to gene symbol mapping (Mandatory).", 
              metavar = "character"),
  
  make_option(c("-o", "--output_dir"), type = "character", default = NULL, 
              help = "Output directory where results will be saved (Mandatory). ", 
              metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# Read inputs

dea = read.delim(opt$input, header = T, sep = "\t")
background_bed = read.delim(opt$background, header = F, sep = "\t")
mapping = read.delim(opt$mapping, header = T, sep = "\t")

dea = merge(dea, mapping, by.x = "gene", by.y = "SYMBOL")

# Get DEGs
for (dea_set in c("all", "upregulated", "downregulated")) {

    if (dea_set == "all") {

        degs = dea %>%
            filter(padj <= 0.05) %>%
            filter(abs(log2FoldChange) >= 2) %>%
            pull(ENTREZID)

    } else if (dea_set == "upregulated") {

        degs = dea %>%
            filter(padj <= 0.05) %>%
            filter(log2FoldChange >= 2) %>%
            pull(ENTREZID)
        
    } else if (dea_set == "downregulated") {

        degs = dea %>%
            filter(padj <= 0.05) %>%
            filter(log2FoldChange <= -2) %>%
            pull(ENTREZID)
        
    }

    out_bed = background_bed %>%
        filter(V4 %in% degs)

    write.table(out_bed, file.path(opt$output_dir, paste0("foreground_", dea_set, ".bed")), col.names = F, row.names = F, quote = F, sep = "\t")

}
