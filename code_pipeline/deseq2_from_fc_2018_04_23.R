# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 4/23/2018
# description: run deseq2 pipeline using FeatureCounts output
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

data_subset = "ribo_w_Female" ## ribo nuclear wc all ribo_w_Female
data_source = "all" ## exon all

## load count matrix  (from count_to_norm_fc.R)
fc = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/fc_gene_from_", 
                    data_source, ".RDS"))

## load previously generated info matrix (fc_count_to_norm.R)
info = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/info.RDS"))

## confirm info row order is same as fc$counts column order
rownames(info) = info$file_name
all(rownames(info) == colnames(fc$counts))


###############################################
# subset of genes detected with all methods
###############################################

# data_subset_list = c("ribo", "nuclear", "wc")
# vsd_f_list = paste0("d1_d2_rnaseq/expression_data_fc/", data_subset_list, "/deseq2_from_", 
#                     data_source, "_vsd2018_08_28.RDS")
# names(vsd_f_list) = data_subset_list
# vsd_list = map(vsd_f_list, ~readRDS(.) %>% rownames)
# gene_subset = Reduce(intersect, vsd_list)

#############################
# normalize counts with 
# DESeq2
#############################

## DESeq2 manual here:
## http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

dds = DESeqDataSetFromMatrix(countData = fc$counts,
                             colData = info,
                             design = ~ Cell_type + gender) # Cell_type + Method + gender

## add length so that fpkm can be calculated
mcols(dds) <- DataFrame(mcols(dds), fc$annotation %>% rename(basepairs = Length, gene = GeneID))

## keeping only a subset based on intron pct or some other metric
# dds = dds[row.names(dds) %in% gene_subset, ]
# ids_to_rm = colnames(dds) %in% c("D1_D2_nuclear.D1.2_L005.sam", "D1_D2_nuclear.D1.5_L005.sam")
# dds = dds[, !ids_to_rm]

## recommended filtering from DESeq2 manual: keep = rowSums(counts(dds)) >= 10
## removing genes with v high counts
low_enough = rowMax(fpkm(dds)) < 5e4 
sum(low_enough)

## removing genes with v low counts (ribo only)
## chosen since minimum lib size is 11.1, so at a minimum 3 fragments per kb 
keep = rowSums(fpkm(dds) > 1/10) >= 2
## for nuclear and WC
keep = rowSums(fpkm(dds) >= 1) >= 1
sum(keep)
keep = keep & low_enough
sum(keep)

dds = dds[keep,]

## running DE
dds = DESeq(dds)
## which contrasts are we running
resultsNames(dds)
saveRDS(dds, paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/deseq2_gene_from_",
                    data_source, "_2018_08_28.RDS"))

# data_subset = "ribo" ## ribo nuclear wc all
# dds = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/deseq2_gene_from_",
#                      data_source, ".RDS"))

## nonstandard chromosomes
# nonstd_chrom_genes = fc$annotation %>% filter(grepl("GL|JH|MT", Chr)) %>% 
#   dplyr::select(GeneID) %>% 
#   filter(GeneID %in% row.names(dds)) %>% 
#   unlist %>% as.character

res = results(dds, name = "Cell_type_D2_vs_D1", independentFiltering = T)
# res = results(dds, name = "Method_wc_vs_nuclear", independentFiltering = T)

summary(res)
## print RPKMs to file (only for genes with non-NA p-values)
non_na_genes = rownames(res)[!is.na(res$padj)]
length(non_na_genes)

fpkm(dds[non_na_genes, ]) %>% as.data.frame %>%
  mutate(gene = non_na_genes) %>%
  select(gene, everything()) %>%
  write_tsv(paste0("d1_d2_rnaseq/de_tables/fc_deseq/rpkm_", data_subset,
                   "_from_", data_source, "_2018_08_28_LENIENT.txt"))

## shrink results
resLFC = lfcShrink(dds, coef="Cell_type_D2_vs_D1") ## specify the name of the coefficient to shrink
resLFC

## look at D1 and D2
drd1 = "ENSMUSG00000021478"
drd2 = "ENSMUSG00000032259"
res[c(drd1, drd2), ]
resLFC[c(drd1, drd2), ]

## summary of results
summary(res)

## plots
# plotMA(res, ylim=c(-2,2)) ## 

## print DEG lists to file
# resOrdered = res[order(res$padj),]
resOrdered = resLFC[order(resLFC$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resSig

resOrdered %>%
  as.data.frame %>% 
  mutate(gene_ens = row.names(resOrdered)) %>%
  filter(!is.na(padj)) %>% 
  select(gene_ens, everything()) %>%
  write_tsv(., paste0("d1_d2_rnaseq/de_tables/fc_deseq/", data_subset,
                      "_d1_v_d2_gene_from_", data_source, "_2018_08_28.txt"))

## plotting counts of single genes
plotCounts(dds, gene=drd1, intgroup="Cell_type") ## Cell_type_D2_vs_D1

## PCA
vsd = varianceStabilizingTransformation(dds, blind=TRUE)
dim(assay(vsd))

saveRDS(assay(vsd), paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/deseq2_from_", 
                           data_source, "_vsd2018_08_28.RDS"))

## for variance partition try both Variance stabilized data and quantlog (same results)
# quantLog = log2( fpm( dds ) + 1)
# saveRDS(quantLog, paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/deseq2_from_", 
#                            data_source, "_quantlog_2018_08_28.RDS"))


### gender Cell_type Method
group_i = "Cell_type"
# DESeq2 is obnoxious and makes you modify the plotPCA function to get PCs beyond 1 and 2
p_data = plotPCA(vsd, intgroup=group_i, returnData=TRUE)
percentVar = round(100 * attr(p_data, "percentVar"), digits = 1)
p = p_data %>%
  ggplot(., aes(PC1, PC2, color=group)) +
  geom_point(size=2.5) + ## 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  # coord_fixed() +
  ### order is: nuclear, ribo, wc
  # scale_color_manual(values = c("royalblue", "red", "black")) + 
  scale_color_manual(values = c("yellow", "blue")) + ##c("blue", "grey")) + ## 
  # geom_text(aes(label = group), col = "black", show.legend = FALSE, check_overlap = F, hjust = "inward") +
  theme_classic()
p

filename = paste0("d1_d2_rnaseq/figures/pca_2018_08_28/", data_subset, 
                  "_from_", data_source, "_pc1_pc2_", group_i, "_2018_08_28.pdf")
ggsave(filename, p, width = 3.5, height = 2.5, units = "in")

#A histogram of pvalues should always be plotted in order to check whether they have been computed correctly.
hist(res$padj, col = "lavender",main = "D1 vs D2", xlab = "p-values")

p_data %>% filter(PC1 > 15)
## D1_D2_ribo.D2F4.sam: abnormal sample, D2 is completely missing and completely off on PCA
## D1_D2_ribo.D1F3.sam: outlier on PCA
p_data %>% filter(PC2 < 2, PC1 < 0)

