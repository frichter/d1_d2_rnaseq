# Felix Richter, Hope Kronman
# Created 6/4/2018
# Description: Get intron/exon %s for all 3
##############################################################

setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("limma", "edgeR", 
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

saveRDS(info, paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/info.RDS"))

saveRDS(fc, paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/fc_gene_from_all.RDS"))

## subset the relevant IDs 
info = readRDS("d1_d2_rnaseq/expression_data_fc/all/info.RDS")
method_list = info$Method %>% unique %>% as.character
data_subset = method_list[[1]] #"nuclear" ## ribo nuclear wc all
subset_ids = info %>% 
  filter(Method == method_list[[1]]) %>% 
  select(file_name) %>% unlist %>% as.character

## load count matrix
fc_all = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/fc_gene_from_all.RDS"))
fc_gene_from_exon = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/fc_gene_from_exon.RDS"))

## confirm filenames are the same for both
identical(fc_all$targets, fc_gene_from_exon$targets)

## create DGEList object
x_all = DGEList(counts = fc_all$counts, genes = fc_all$annotation)
x_gene_from_exon = DGEList(counts = fc_gene_from_exon$counts, genes = fc_gene_from_exon$annotation)

## filter out low-expression features (threshold may need adjustment)
## relaxed: keep if RNA cpm is greater than 1 in at least 2 samples (both abritrary)
isexpr = rowSums(cpm(x_all) > 1) >= 2
## stringent: keep if avg RPKM across all samples is >1
# isexpr = rowMeans(rpkm(x)) >= 1
x_all = x_all[isexpr,]
x_gene_from_exon = x_gene_from_exon[isexpr,]

# confirm everything is positive
counts_introns = x_all$counts - x_gene_from_exon$counts
any(counts_introns < 0)
all(counts_introns >= 0)

# percent introns/exons:
pct_exon = x_gene_from_exon$counts/x_all$counts
pct_intron = 1 - pct_exon
dim(pct_exon)
pct_exon[1:5, 1:5]
pct_intron[1:5, 1:5]

mean_pct_exon = rowMeans(pct_exon[, subset_ids], na.rm = T)
sd_pct_exon = apply(pct_exon[, subset_ids], 1, sd, na.rm = T)

pct_exon_df = as.data.frame(mean_pct_exon) %>% 
  mutate(sd_pct_exon = sd_pct_exon,
         ci_hi = mean_pct_exon + sd_pct_exon*1.96,
         ci_lo = mean_pct_exon - sd_pct_exon*1.96) %>% 
  mutate(ci_lo = ifelse(ci_lo < 0, 0, ci_lo)) %>% 
  mutate(gene = row.names(pct_exon)) %>% 
  mutate(Method = data_subset)

p = ggplot(pct_exon_df, aes(x = mean_pct_exon, y = mean_pct_exon)) +
  # geom_density() + 
  geom_point() + 
  # possible interesting options: geom_freqpoly, geom_rug, geom_path()
  # geom_density_2d and/or stat_density_2d(aes(fill = ..level..), geom="polygon")
  # geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.5, fill = "grey50") +
  theme_classic()
p



