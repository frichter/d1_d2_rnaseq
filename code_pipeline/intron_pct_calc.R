# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 4/23/2018
# description: compare exon and intron percents
##############################################################

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("limma", "edgeR", "annotate", "org.Mm.eg.db", "DESeq2",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

#############################
# Load data
#############################

data_subset = "all" ## ribo nuclear wc all ribo_w_Female
data_source = "all" ## exon all

## load count matrix  (from count_to_norm_fc.R)
fc_all = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/fc_gene_from_", 
                        data_source, ".RDS"))
# fc_gene_from_all.RDS fc_gene_from_exon.RDS
data_source = "exon" ## exon all
fc_exon = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/fc_gene_from_", 
                        data_source, ".RDS"))

fc_all %>% names
fc_all$counts %>% dim
fc_exon$counts %>% dim

# fc_intron = fc_all$counts - fc_exon$counts
# ## confirm none are less than 0
# sum(fc_intron < 0)

## load previously generated info matrix
info = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/info.RDS"))

## confirm info row order is same as fc$counts column order
rownames(info) = info$file_name
all(rownames(info) == colnames(fc_all$counts))

## get lists of genes in the 3 methods
method_vec = c("nuclear", "wc", "ribo")
names(method_vec) = method_vec
kept_genes = map(method_vec, ~ readRDS(
  paste0("d1_d2_rnaseq/expression_data_fc/", .,
         "/deseq2_from_all_vsd2018_07_12.RDS")) %>% row.names)
# common_genes = Reduce(intersect, kept_genes)
kept_genes_df = map2_df(kept_genes, names(kept_genes),
                     function(x, y) cbind("ens_gene" = x, "Method" = y) %>% as.data.frame)
kept_genes_df %>% group_by(Method) %>% tally

###########################
# combine exon and all data
###########################

all_long = fc_all$counts %>% as.data.frame %>% 
  mutate(ens_gene = row.names(fc_all$counts )) %>% 
  gather(key = "file_name", value = "ct_all", -ens_gene)

exon_long = fc_exon$counts %>% as.data.frame %>% 
  mutate(ens_gene = row.names(fc_exon$counts )) %>% 
  gather(key = "file_name", value = "ct_exon", -ens_gene)

cts_long = inner_join(all_long, exon_long)
cts_long %<>% inner_join(info %>% select(file_name, Method)) %>% inner_join(kept_genes_df)

###########################
# calculating % intronic
###########################

pct_intron = cts_long %>% 
  mutate(ct_intron = ct_all - ct_exon) %>%
  # mutate(pct_intron = ((ct_all - ct_exon)/ct_all)*100) %>% 
  group_by(Method, ens_gene) %>% 
  summarise(ct_all = mean(ct_all), ct_exon = mean(ct_exon), ct_intron = mean(ct_intron)) %>% 
  ungroup %>% 
  mutate(pct_intron = 100 * (ct_intron/ct_all))

pct_intron %>% 
  select(Method, ens_gene, ct_all, pct_intron) %>% 
  write_tsv("d1_d2_rnaseq/de_tables/det_inconsistency_exploration/pct_intron.txt")

p = pct_intron %>% 
  mutate(gene_interest = ens_gene %in% overlap_transcripts) %>% 
  group_by(Method, gene_interest) %>% mutate(med = median(pct_intron)) %>% ungroup %>% 
  mutate(Method = factor(Method, levels = c("nuclear", "wc", "ribo"))) %>% 
  ggplot(aes(x = Method, y = pct_intron, color = gene_interest)) +
  geom_violin(adjust = 0.5, scale = "width", trim = T) +
  geom_boxplot(aes(y = med), size = 1) +
  # scale_color_manual(values = c("royalblue", "red", "black")) +
  theme_classic() +
  xlab("") + ylab("Percent intronic reads\n(average per gene)")
p
ggsave("d1_d2_rnaseq/figures/intronic_reads_pct/intronic_read_pct_2018_08_24.png",
       p, width = 3.2, height = 2.8)

###########################
# just looking at DRD1/DRD2
###########################

drd1 = "ENSMUSG00000021478"
drd2 = "ENSMUSG00000032259"

pct_drd = cts_long %>% 
  filter(ens_gene %in% c(drd1, drd2)) %>% 
  mutate(ct_intron = ct_all - ct_exon) %>%
  mutate(pct_intron = 100 * (ct_intron/ct_all)) %>% 
  group_by(Method, ens_gene) %>% mutate(mid = median(pct_intron)) %>% ungroup %>% 
  mutate(gene = ifelse(ens_gene == drd1, "Drd1", "Drd2"))

p = pct_drd %>% 
  mutate(Method = factor(Method, levels = c("nuclear", "wc", "ribo"))) %>% 
  ggplot(aes(x = Method, y = pct_intron)) +
  # geom_violin() +
  geom_point(size = 0.75) +
  geom_boxplot(aes(y = mid, col = Method), size = 1) +
  scale_color_manual(values = c("royalblue", "red", "black")) +
  facet_wrap(~gene) +
  theme_classic() +
  ylim(0, 100) +
  xlab("") + ylab("Percent intronic reads")
p
ggsave("d1_d2_rnaseq/figures/intronic_reads_pct/drd_intronic_read_pct_2018_08_24.png",
       p, width = 4, height = 2.5)


