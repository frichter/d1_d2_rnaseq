# Felix Richter, Hope Kronman
# 8/10/2017
# normalize count matrix RNAseq output, and integrate metadata

## set the home directory
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

#############################
# FeatureCounts data cleaning
#############################

## load count matrices
cts_file_list = list.files("d1_d2_rnaseq/expression_data_fc/counts_matrices", "counts_.*.RDS", full.names = T)
names(cts_file_list) = gsub(".*matrices/counts_", "", cts_file_list) %>% gsub("_17_11.RDS", "", .)
cts_list = map(cts_file_list, readRDS)

fc = readRDS("d1_d2_rnaseq/expression_data_fc/counts_matrices/counts_all_17_11.RDS")

## print file names to confirm you are looking at the correct data
fc$targets

## create metadata matrix from filenames 
all_targets = fc$targets
info = as.data.frame(all_targets) %>% 
  rename(file_name = all_targets) %>% 
  mutate(sub_name = gsub("D1_D2_", "", file_name) %>% gsub(".sam", "", .) %>% 
           gsub("CTRL", "_", .) %>% gsub("whole_cell", "wc", .) %>% 
           gsub("_repeat", "Rpt", .)) %>% 
  separate(sub_name, c("Method", "Cell_type", "Replicate", "Lane"), fill = "right")

## how many samples are there for each?
info %>% group_by(Method, Cell_type) %>% tally

#############################
# plot alignment statistics: 
# are there any outliers/abnormal samples?
#############################

stats_long = fc$stat %>% as.data.frame %>% 
  gather("file_name", "read_ct", -Status) %>% 
  ## remove lines without reads
  filter(read_ct != 0) %>% 
  ## calculate percentages
  group_by(file_name) %>% 
  mutate(read_pct = 100 * read_ct/sum(read_ct)) %>% 
  ungroup %>%  
  ## join with info
  left_join(info) 

p = stats_long %>% 
  ggplot(., aes(x = read_ct, fill = Method)) + 
  geom_histogram(bins = 25) +
  facet_wrap(~ Status) +
  scale_fill_manual(values = c("purple3", "skyblue", "darkblue")) + 
  theme_classic()
p
ggsave("d1_d2_rnaseq/figures/qc_fc/mapping_features_ct.png", p, width = 5, height = 2.5)

## are there any counts from same sample but different lanes? If yes then could add the counts
filenames_to_add = info %>% group_by(Method, Cell_type, Replicate) %>%
  mutate(n = n()) %>% ungroup %>% filter(n > 1) %>% select(filename) %>% unlist %>% as.character
fc$counts[, filenames_to_add] %>% head
# new_name = gsub("L00[5-8]", "Lx", filenames_to_add) %>% unique
# new_cts = rowSums(fc$counts[, filenames_to_add])
# fc$counts %>% head
# fc$counts %<>% cbind(., new_cts)
# colnames(fc$counts)[16] = new_name


## prepare info matrix
## clean factor levels
info %<>% mutate_all(as.factor) %>% mutate_all(droplevels)
saveRDS(info, "d1_d2_rnaseq/expression_data_fc/all_info.RDS")

## create DGEList object
x = DGEList(counts = fc$counts, genes = fc$annotation)

## filter out low-expression features (threshold may need adjustment)
## relaxed: keep if RNA cpm is greater than 1 in at least 2 samples (both abritrary)
isexpr = rowSums(cpm(x) > 1) >= 2
## stringent: keep if avg RPKM across all samples is >1
# isexpr = rowMeans(rpkm(x)) >= 1
x = x[isexpr,]

## normalize to library size within and between experiments
x = calcNormFactors(x)

## save count data
saveRDS(x, "d1_d2_rnaseq/expression_data_fc/all_norm_strict.RDS") 
## nuc_norm.RDS wc_norm.RDS all_norm.RDS allD1_norm.RDS allD2_norm.RDS


#########################
# HTSeq data cleaning
#########################

nuclear_cts = read_tsv("d1_d2_rnaseq/expression_data_htseq/counts/NUCLEAR_Baseline_counts_mm10.txt")
whole_cell_cts = read_tsv("d1_d2_rnaseq/expression_data_htseq/counts/WHOLE_CELL_Baseline_counts_M13.txt")

## clean column names
names(nuclear_cts) %<>% gsub("Sample_", "", .) %>% gsub("\\.htseq_counts.txt", "", .) %>% make.names

## prepare metadata matrices
info_nuc = PrepNucMetadataHTS(nuclear_cts) 
saveRDS(info_nuc, "d1_d2_rnaseq/expression_data_htseq/nuc_info.RDS")
info_wc = PrepWCMetadataHTS(whole_cell_cts)
saveRDS(info_wc, "d1_d2_rnaseq/expression_data_htseq/wc_info.RDS")

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
