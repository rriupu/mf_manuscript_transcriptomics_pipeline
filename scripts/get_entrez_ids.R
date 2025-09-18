library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(optparse)

## -------------- ##
## Read arguments ##
## -------------- ##

option_list = list(

  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "BED file with the raw TSS annotations (Mandatory).", 
              metavar = "character"),
  
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "Path to the output BED file (Mandatory). ", 
              metavar = "character")
)

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# Read input files

bed = read.delim(opt$input, header = F, sep = "\t")
bed = bed %>%
    filter(str_detect(V4, "NM_")) # NM_ ids correspond to mRNAs

refseq_ids = sapply(strsplit(bed$V4, ".", fixed = T), "[[", 1)
bed$V4 = refseq_ids

# Annotations
entrez = AnnotationDbi::select(org.Hs.eg.db, keys = refseq_ids, columns = c("REFSEQ", "ENTREZID"), keytype = "REFSEQ")

bed = merge(bed, entrez, by.x = "V4", by.y = "REFSEQ") %>%
    filter(ENTREZID != "")
bed = bed[, c("V1", "V2", "V3", "ENTREZID", "V5", "V6")]
bed = bed[!is.na(bed$ENTREZID), ]

write.table(bed, file = opt$output, quote = FALSE, col.names = F, row.names = F, sep = "\t")