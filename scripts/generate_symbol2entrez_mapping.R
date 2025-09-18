library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(optparse)

option_list = list(

  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "TSV file with one of the DEA results (Mandatory).", 
              metavar = "character"),
  
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "Output mapping (Mandatory). ", 
              metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

dea = read.delim(opt$input, header = T, sep = "\t")
genes = dea$gene

entrez = AnnotationDbi::select(org.Hs.eg.db, keys = genes, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
entrez = entrez[!is.na(entrez$ENTREZID), ]

write.table(entrez, opt$output, col.names = T, row.names = F, quote = F, sep = "\t")