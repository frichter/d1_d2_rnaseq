# Felix Richter, Hope Kronman
# 8/10/2017
# normalize count matrix RNAseq output, and integrate metadata

setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("limma", "edgeR", "annotate", "org.Mm.eg.db", 
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

## load custom functions
source("d1_d2_rnaseq/code_pipeline/count_to_norm_functions.R")

## load count matrices
nuclear_cts = read_tsv("d1_d2_rnaseq/expression_data/counts/NUCLEAR_Baseline_counts_mm10.txt")
whole_cell_cts = read_tsv("d1_d2_rnaseq/expression_data/counts/WHOLE_CELL_Baseline_counts_M13.txt")

## clean column names
names(nuclear_cts) %<>% gsub("Sample_", "", .) %>% gsub("\\.htseq_counts.txt", "", .) %>% make.names

## prepare metadata matrices
info_nuc = PrepNucMetadata(nuclear_cts) 
saveRDS(info_nuc, "d1_d2_rnaseq/expression_data/nuc_info.RDS")
info_wc = PrepWCMetadata(whole_cell_cts)
saveRDS(info_wc, "d1_d2_rnaseq/expression_data/wc_info.RDS")

info_all = bind_rows("Nuc" = info_nuc, "WC" = info_wc, .id = "Method")
## clean factor levels
info_all %<>% mutate_all(as.factor) %>% mutate_all(droplevels)
saveRDS(info_all, "d1_d2_rnaseq/expression_data/all_info.RDS")

## create combined matrices
all_cts = inner_join(nuclear_cts, whole_cell_cts, by = "gene_id")
## create D1 matrix
d1_ids = info_all %>% filter(Receptor == "D1") %>% select(ID) %>% unlist %>% as.character
allD1_cts = all_cts[, c("gene_id", d1_ids)]
info_all %>% filter(Receptor == "D1") %>% saveRDS("d1_d2_rnaseq/expression_data/allD1_info.RDS")
## create D2 matrix
d2_ids = info_all %>% filter(Receptor == "D2") %>% select(ID) %>% unlist %>% as.character
allD2_cts = all_cts[, c("gene_id", d2_ids)]
info_all %>% filter(Receptor == "D2") %>% saveRDS("d1_d2_rnaseq/expression_data/allD2_info.RDS")

## convert to matrix
ct_mtx = allD2_cts %>% select(-gene_id) %>% as.matrix ## nuclear_cts whole_cell_cts all_cts
row.names(ct_mtx) = allD2_cts$gene_id ## nuclear_cts whole_cell_cts all_cts

## create DGEList object
x <- DGEList(counts = ct_mtx)

## filter out low-expression features (threshold may need adjustment)
## keep if RNA cpm is greater than 1 in at least 2 samples (both abritrary)
isexpr <- rowSums(cpm(x) > 1) >= 2
x <- x[isexpr,]

## normalize to library size within and between experiments
x <- calcNormFactors(x)

## save count data
saveRDS(x, "d1_d2_rnaseq/expression_data/allD2_norm.RDS") 
## nuc_norm.RDS wc_norm.RDS all_norm.RDS allD1_norm.RDS allD2_norm.RDS


######################
# One time use
######################

## check overlap in gene IDs between nuclear_cts and whole_cell_cts
nuclear_cts %>% filter(!(gene_id %in% whole_cell_cts$gene_id)) %>% dim
whole_cell_cts %>% filter(!(gene_id %in% nuclear_cts$gene_id)) %>% dim

all_ensIDs = c(whole_cell_cts$gene_id, nuclear_cts$gene_id) %>% unique

GetGeneInfo = function(ensIDs) {
  # add gene symbol information
  geneSyms = AnnotationDbi::select(org.Mm.eg.db, ensIDs, "SYMBOL","ENSEMBL")
  # if gene symbol not present, keep ID with Unk: prefix
  geneSyms %<>% mutate(SYMBOL = ifelse(is.na(SYMBOL), paste0("Unk:",ENSEMBL), SYMBOL))
  return(geneSyms)
}
## import gene annotations

gene_df = GetGeneInfo(row.names(ct_mtx))
gene_df %>% filter(duplicated(ENSEMBL))
