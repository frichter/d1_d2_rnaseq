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

data_subset = "ribo" ## ribo nuclear wc all

## load count matrix  (from count_to_norm_fc.R)
fc = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/fc_gene_from_exon.RDS"))
# fc_gene_from_all.RDS

## load previously generated info matrix
info = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/info.RDS"))

## confirm info row order is same as fc$counts column order
rownames(info) = info$file_name
all(rownames(info) == colnames(fc$counts))

#############################
# normalize counts with 
# DESeq2
#############################

## DESeq2 manual here:
## http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

dds = DESeqDataSetFromMatrix(countData = fc$counts,
                             colData = info,
                             design = ~ Cell_type) # Cell_type + Method + gender
dds

## recommended filtering from DESeq2 manual: keep = rowSums(counts(dds)) >= 10
## filtering based on fragments per million
keep = rowSums(fpm(dds) > 1) >= 2 ## DESeq2 also has fpkm
sum(keep)
# how many genes are we keeping, ~ 30k for nuclear, 34k for wc, 12k for ribo
# For counts from all (not just from exons): wc 34510, nuclear 28476, ribo 11133, all 28476
# For counts from exons: wc 39629, nuclear 38387, ribo 10545, all 
dds = dds[keep,]

## running DE
dds = DESeq(dds)
## which contrasts are we running
resultsNames(dds)
# saveRDS(dds, "d1_d2_rnaseq/expression_data_fc/all/deseq2_gene_from_all.RDS")
saveRDS(dds, paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/deseq2_gene_from_exon.RDS"))

res = results(dds, name = "Cell_type_D2_vs_D1")
res

## shrink results
resLFC = lfcShrink(dds, coef=2)
resLFC

## look at D1 and D2
drd1 = "ENSMUSG00000021478"
drd2 = "ENSMUSG00000032259"
res[c(drd1, drd2), ]
resLFC[c(drd1, drd2), ]

## adjusted p-values were less than 0.05. 
## 14k genes that are DE for nuclear?? Only 612 for WC (which makes more sense..)
## 1k for riboseq
sum(resLFC$padj < 0.05, na.rm=TRUE)
summary(res)

## plots
plotMA(resLFC, ylim=c(-2,2)) ## 

# identify outliers
res05 <- results(dds, alpha=0.05)
summary(res05)

## print DEG lists to file
resOrdered = res[order(resLFC$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resSig

resSig %>% as.data.frame %>% 
  mutate(gene_ens = row.names(resSig)) %>% 
  select(gene_ens, everything()) %>% 
  write_tsv(., paste0("d1_d2_rnaseq/de_tables/fc_deseq/", data_subset, "_d1_v_d2_gene_from_exon_2018_05_05.txt"))
# _d1_v_d2_gene_from_all_2018_05_05.txt _d1_v_d2_gene_from_exon_2018_05_05.txt

## plotting counts of single genes
plotCounts(dds, gene=drd1, intgroup="Cell_type") ## Cell_type_D2_vs_D1

##
vsd = vst(dds, blind=FALSE)
p_data = plotPCA(vsd, intgroup="Cell_type", returnData=TRUE)
percentVar = round(100 * attr(p_data, "percentVar"), digits = 1)
p = p_data %>% 
  left_join(info %>% select(file_name, ID, gender, Method), by = c("name" = "file_name")) %>% 
  ggplot(., aes(PC1, PC2, color=Cell_type)) + ## gender Cell_type Method
  geom_point() + ## size=3
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  # coord_fixed() +
  # scale_color_manual(values = c("skyblue", "red3", "purple3")) + 
  scale_color_manual(values =  c("green3", "red")) + ##c("blue", "grey")) + ##
  # geom_text(aes(label = ID), col = "black", show.legend = FALSE, check_overlap = F, hjust = "inward") + 
  theme_classic()
p
filename = paste0("d1_d2_rnaseq/figures/qc_fc/deseq2/", data_subset, "_pc1_pc2_gender_2018_01_29.png")
## _pc1_pc2_gender.png _pc1_pc2_cell_type.png
ggsave(filename, p, width = 4, height = 3, units = "in")

#A histogram of p–values should always be plotted in order to check whether they have been computed correctly.
hist(res$padj, col = "lavender",main = "D1 vs D2", xlab = "p-values")

p_data %>% filter(PC1 > 15)
## D1_D2_ribo.D2F4.sam: abnormal sample, D2 is completely missing and completely off on PCA
## D1_D2_ribo.D1F3.sam: outlier on PCA


########################################################
# obtain variance stabilized expression for WGCNA
########################################################

## DESeq2 manual here:
## http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## confirm info row order is same as fc$counts column order
rownames(info) = info$file_name
all(rownames(info) == colnames(fc$counts))

dds = DESeqDataSetFromMatrix(countData = fc$counts,
                             colData = info,
                             design = ~ Cell_type + gender) # Cell_type + Method + gender
dds

## recommended filtering from DESeq2 manual: keep = rowSums(counts(dds)) >= 10
## filtering based on fragments per million
keep = rowSums(fpm(dds) > 1) >= 2 ## DESeq2 also has fpkm
sum(keep) ## how many genes are we keeping, ~ 30k for nuclear, 34k for wc, 12k for ribo
dds = dds[keep,]
# vsd = vst(dds, blind=TRUE)

vsd = varianceStabilizingTransformation(dds, blind=TRUE)
saveRDS(assay(vsd), paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/deseq2_vsd2018_04_23.RDS"))

# vsd = varianceStabilizingTransformation(dds, blind=TRUE, fitType="parametric")
# vst is a wrapper around varianceStabilizingTransformation where  blind=FALSE
# vsd = vst(dds, blind=FALSE)
# blind estimation of variance function is recommended according to docs.
# Not sure if blind=F would regress out design matrix effects
## possibly collapseReplicates to combine counts from technical replicates


