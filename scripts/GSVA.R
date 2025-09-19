## --------------- ##
## Library loading ##
## --------------- ##

library(SummarizedExperiment)
library(optparse)
library(edgeR)
library(org.Hs.eg.db)
library(GO.db)
library(GSVA)
library(limma)
library(sva)
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

## ------------ ##
## Prepare data ##
## ------------ ##

dset = prepareDESeqDataset(conditionMapping, geneCountsDir, formula = "~ donor + condition")

smallestGroupSize = 4
keep = rowSums(counts(dset) >= 10) >= smallestGroupSize
dset = dset[keep, ]

dds = DESeq(dset)

# Get CPMs
lcpms = cpm(counts(dds, normalized = T), log = TRUE)
coldata = colData(dds)

se = SummarizedExperiment(assays = list(counts = lcpms), colData = coldata)

## ------------------- ##
## Prepare annotations ##
## ------------------- ##

goannot = AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns=c("GO", "SYMBOL", "ONTOLOGY"))
goannot = goannot[goannot$ONTOLOGY == "BP", ]
go_names = AnnotationDbi::select(GO.db, keys = goannot$GO, columns = c("TERM", "GOID"))

genesbygo = split(goannot$SYMBOL, goannot$GO)

## -------- ##
## Run GSVA ##
## -------- ##

pars = gsvaParam(se, genesbygo, minSize = 5, maxSize = 500)
res = gsva(pars)

mod = model.matrix(~ condition, colData(res))
mod0 = model.matrix(~ 1, colData(res))
sv = sva(assay(res), mod, mod0)
mod = cbind(mod, sv$sv)
fit = lmFit(assay(res), mod)
fit.eb = eBayes(fit, robust = TRUE)
dea = decideTests(fit.eb)
summary(dea)

gssizes = geneSetSizes(res)
plot(sqrt(gssizes), sqrt(fit.eb$sigma), xlab = "Sqrt(gene sets sizes)",
          ylab = "Sqrt(standard deviation)", las = 1, pch = ".", cex = 4)
lines(lowess(sqrt(gssizes), sqrt(fit.eb$sigma)), col = "red", lwd = 2)

fit.eb.trend = eBayes(fit, robust = TRUE, trend = gssizes)
dea = decideTests(fit.eb.trend)
summary(dea)

## --------------------- ##
## Get enriched pathways ##
## --------------------- ##

tt = topTable(fit.eb.trend, coef = 2, n = Inf)
DEpwys = rownames(tt)[tt$adj.P.Val <= 0.05]
length(DEpwys)
names(DEpwys) = sapply(
    DEpwys,
    function(x) {
        term = go_names %>%
            filter(GOID == x) %>%
            distinct() %>%
            pull(TERM)
        return(term)
    }
)

DEpwys_es = removeBatchEffect(
    x = assay(res[DEpwys, ]),
    covariates = mod[, 24:ncol(mod)],
    design = mod[, 1:23]
)
rownames(DEpwys_es) = names(DEpwys)

heatmap_mat = assay(res[DEpwys, ])

rownames(heatmap_mat) = names(DEpwys)

## ----------------------- ##
## Focus on relevant terms ##
## ----------------------- ##

terms = c(
    "positive regulation of chemokine (C-X-C motif) ligand 2 production",
    "positive regulation of T cell activation",
    "positive chemotaxis",
    "T-helper 1 type immune response",
    "myeloid cell apoptotic process",
    "positive regulation of potassium ion transmembrane transport",
    "peptide antigen assembly with MHC class II protein complex",
    "positive regulation of immune response",
    "antigen processing and presentation of exogenous peptide antigen via MHC class II",
    "positive regulation of CD4-positive, CD25-positive, alpha-beta regulatory T cell differentiation",
    "positive regulation of memory T cell differentiation",
    "positive regulation of MHC class II biosynthetic process",
    "positive regulation of cytokine-mediated signaling pathway",
    "positive regulation of fibroblast migration",
    "negative regulation of SMAD protein signal transduction",
    "cell migration",
    "blood vessel remodeling",
    "cell chemotaxis",
    "positive regulation of T cell migration",
    "negative regulation of type 2 immune response",
    "positive regulation of interleukin-23 production",
    "negative regulation of T cell apoptotic process",
    "positive regulation of T cell chemotaxis",
    "positive regulation of neutrophil apoptotic process",
    "granulocyte chemotaxis",
    "wound healing, spreading of cells",
    "negative regulation of T-helper 2 cell differentiation",
    "establishment of T cell polarity",
    "positive regulation of inflammatory response",
    "positive regulation of non-canonical NF-kappaB signal transduction",
    "defense response to bacterium",
    "phagocytosis",
    "chemotaxis",
    "positive regulation of angiogenesis",
    "positive regulation of BMP signaling pathway",
    "immune response-regulating signaling pathway",
    "positive regulation of CD4-positive, alpha-beta T cell differentiation",
    "natural killer cell degranulation",
    "macrophage activation",
    "interleukin-10-mediated signaling pathway",
    "negative regulation of T-helper 17 cell differentiation"
)

## ----------------- ##
## Generate heatmaps ##
## ----------------- ##

heatmap_mat = heatmap_mat[rownames(heatmap_mat) %in% terms, ]

heatmap = Heatmap(
    heatmap_mat,
    column_km = 9,
    row_km = 3,
    top_annotation = HeatmapAnnotation(
        Condition = res$condition,
        Donor = res$donor
    ),
    show_row_names = T,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 8)
)

pdf(file.path(outputDir, "GO_heatmap_mat.pdf"), width = 14, height = 10)
draw(heatmap)
dev.off()

h = Heatmap(
    DEpwys_es[rownames(DEpwys_es) %in% terms, ],
    column_km = 9,
    row_km = 3,
    top_annotation = HeatmapAnnotation(
        Condition = res$condition,
        Donor = res$donor
    ),
    show_row_names = T,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 8)
)
pdf(file.path(outputDir, "GSVA_GO.pdf"), width = 14, height = 10)
draw(h)
dev.off()

## ------------------- ##
## Reactome enrichment ##
## ------------------- ##

# URL = "https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Hs/c2.cp.reactome.v2025.1.Hs.symbols.gmt"

c2_reactome = clusterProfiler::read.gmt("scripts/c2.cp.reactome.v2025.1.Hs.symbols.gmt")
terms = unique(c2_reactome$term)
genestopath = lapply(terms, function(x) {
    genes = unique(c2_reactome$gene[c2_reactome$term == x])
    return(genes)
})
names(genestopath) = terms

pars = gsvaParam(se, genestopath, minSize = 5, maxSize = 500)
res = gsva(pars)


mod = model.matrix(~ condition, colData(res))
mod0 = model.matrix(~ 1, colData(res))
sv = sva(assay(res), mod, mod0)
mod = cbind(mod, sv$sv)
fit = lmFit(assay(res), mod)
fit.eb = eBayes(fit, robust = TRUE)
dea = decideTests(fit.eb)
summary(dea)

gssizes = geneSetSizes(res)
plot(sqrt(gssizes), sqrt(fit.eb$sigma), xlab = "Sqrt(gene sets sizes)",
          ylab = "Sqrt(standard deviation)", las = 1, pch = ".", cex = 4)
lines(lowess(sqrt(gssizes), sqrt(fit.eb$sigma)), col = "red", lwd = 2)

fit.eb.trend = eBayes(fit, robust = TRUE, trend = gssizes)
dea = decideTests(fit.eb.trend)
summary(dea)

tt = topTable(fit.eb.trend, coef = 2, n = Inf)
DEpwys = rownames(tt)[tt$adj.P.Val <= 0.05]
length(DEpwys)

DEpwys_es = removeBatchEffect(
    x = assay(res[DEpwys, ]),
    covariates = mod[, 24:ncol(mod)],
    design = mod[, 1:23]
)

terms = c(
    "REACTOME_INFLUENZA_INFECTION",
    "REACTOME_SARS_COV_1_HOST_INTERACTIONS",
    "REACTOME_SARS_COV_1_MODULATES_HOST_TRANSLATION_MACHINERY",
    "REACTOME_MTOR_SIGNALLING",
    "REACTOME_SIGNALING_BY_FGFR2",
    "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES",
    "REACTOME_SIGNALING_BY_FGFR2_IN_DISEASE",
    "REACTOME_SODIUM_CALCIUM_EXCHANGERS",
    "REACTOME_SIGNALING_BY_FGFR2_IIIA_TM",
    "REACTOME_SIGNALING_BY_GPCR",
    "REACTOME_REGULATION_OF_TP53_EXPRESSION_AND_DEGRADATION",
    "REACTOME_GPCR_LIGAND_BINDING",
    "REACTOME_SARS_COV_1_INFECTION",
    "REACTOME_CD163_MEDIATING_AN_ANTI_INFLAMMATORY_RESPONSE",
    "REACTOME_SARS_COV_2_INFECTION",
    "REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS",
    "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_G2_CELL_CYCLE_ARREST",
    "REACTOME_REGULATION_OF_TP53_ACTIVITY",
    "REACTOME_APC_CDC20_MEDIATED_DEGRADATION_OF_NEK2A",
    "REACTOME_DNA_DAMAGE_RECOGNITION_IN_GG_NER",
    "REACTOME_GAB1_SIGNALOSOME",
    "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DNA_REPAIR_GENES",
    "REACTOME_G0_AND_EARLY_G1",
    "REACTOME_VIRAL_MESSENGER_RNA_SYNTHESIS",
    "REACTOME_AMINO_ACIDS_REGULATE_MTORC1",
    "REACTOME_SARS_COV_INFECTIONS",
    "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING",
    "REACTOME_RRNA_MODIFICATION_IN_THE_MITOCHONDRION",
    "REACTOME_INTERLEUKIN_10_SIGNALING",
    "REACTOME_APC_C_CDC20_MEDIATED_DEGRADATION_OF_CYCLIN_B",
    "REACTOME_INTERLEUKIN_20_FAMILY_SIGNALING",
    "REACTOME_GROWTH_HORMONE_RECEPTOR_SIGNALING",
    "REACTOME_MTORC1_MEDIATED_SIGNALLING",
    "REACTOME_G1_S_SPECIFIC_TRANSCRIPTION",
    "REACTOME_INHIBITION_OF_THE_PROTEOLYTIC_ACTIVITY_OF_APC_C_REQUIRED_FOR_THE_ONSET_OF_ANAPHASE_BY_MITOTIC_SPINDLE_CHECKPOINT_COMPONENTS",
    "REACTOME_REGULATION_OF_TP53_ACTIVITY_THROUGH_PHOSPHORYLATION",
    "REACTOME_ABERRANT_REGULATION_OF_MITOTIC_EXIT_IN_CANCER_DUE_TO_RB1_DEFECTS",
    "REACTOME_SIGNALING_BY_FGFR",
    "REACTOME_REGULATION_OF_INSULIN_LIKE_GROWTH_FACTOR_IGF_TRANSPORT_AND_UPTAKE_BY_INSULIN_LIKE_GROWTH_FACTOR_BINDING_PROTEINS_IGFBPS",
    "REACTOME_OTHER_INTERLEUKIN_SIGNALING",
    "REACTOME_DISEASES_OF_MITOTIC_CELL_CYCLE",
    "REACTOME_REGULATION_OF_IFNA_IFNB_SIGNALING",
    "REACTOME_INTERLEUKIN_21_SIGNALING",
    "REACTOME_INTERLEUKIN_9_SIGNALING",
    "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_CYCLE_GENES",
    "REACTOME_REGULATION_OF_TP53_ACTIVITY_THROUGH_METHYLATION",
    "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS",
    "REACTOME_S_PHASE",
    "REACTOME_SIGNALING_BY_INTERLEUKINS",
    "REACTOME_EXTENSION_OF_TELOMERES",
    "REACTOME_MITOTIC_SPINDLE_CHECKPOINT",
    "REACTOME_REGULATION_OF_IFNG_SIGNALING",
    "REACTOME_EPITHELIAL_MESENCHYMAL_TRANSITION_EMT_DURING_GASTRULATION",
    "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
    "REACTOME_IFNG_SIGNALING_ACTIVATES_MAPKS",
    "REACTOME_ABERRANT_REGULATION_OF_MITOTIC_G1_S_TRANSITION_IN_CANCER_DUE_TO_RB1_DEFECTS",
    "REACTOME_WNT_MEDIATED_ACTIVATION_OF_DVL",
    "REACTOME_SIGNALING_BY_LTK_IN_CANCER"
)

heatmap_mat = DEpwys_es[rownames(DEpwys_es) %in% terms, ]
rownames(heatmap_mat) = gsub("REACTOME_", "", rownames(heatmap_mat))

h = Heatmap(
    heatmap_mat,
    column_km = 9,
    row_km = 3,
    top_annotation = HeatmapAnnotation(
        Condition = res$condition,
        Donor = res$donor
    ),
    show_row_names = T,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 8)
)
pdf(file.path(outputDir, "GSVA_reactome.pdf"), width = 14, height = 10)
draw(h)
dev.off()

## ------------------------------------ ##
## Metapatient as average per condition ##
## ------------------------------------ ##

conditions_ls = list()
for (i in seq_along(unique(conditionMapping$condition))) {

    condition = unique(conditionMapping$condition)[i]
    conditions_ls[[i]] = as.data.frame(rowMeans(heatmap_mat[, grepl(paste0("_", gsub("_", "-", condition), "_", sep = ""), colnames(heatmap_mat))]))
    colnames(conditions_ls[[i]]) = condition
    names(conditions_ls)[i] = condition

}
df = do.call(cbind, conditions_ls)
h2 = Heatmap(
    as.matrix(df),
    column_km = 2,
    row_km = 2,
    show_row_names = T,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 8)
)
pdf(file.path(outputDir, "GSVA_metapatient.pdf"), width = 14, height = 10)
draw(h2)
dev.off()
