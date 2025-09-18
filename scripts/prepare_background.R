library(dplyr)
library(optparse)

## -------------- ##
## Read arguments ##
## -------------- ##

option_list = list(

  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "TSV file with one of the DEA results (Mandatory).", 
              metavar = "character"),

  make_option(c("-b", "--bed"), type = "character", default = NULL, 
              help = "BED file with all the promoter coordinates (Mandatory).", 
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
bed = read.delim(opt$bed, header = F, sep = "\t")
mapping = read.delim(opt$mapping, header = T, sep = "\t")

dea = merge(dea, mapping, by.x = "gene", by.y = "SYMBOL")

out_bed = bed %>%
    filter(V4 %in% unique(dea$ENTREZID))

write.table(out_bed, file.path(opt$output_dir, "background.bed"), col.names = F, row.names = F, quote = F, sep = "\t")
