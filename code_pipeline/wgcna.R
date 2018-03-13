# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 3/13/2018
# description: WGCNA
##############################################################

# export TMP=/sc/orga/projects/chdiTrios/Felix/dna_rna/eqtl_wgs/fastqtl_2018_01/
# module load R/3.4.1 ## R/3.3.1 ## use this minerva version for for limma
# R

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("WGCNA", 
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr", "colorout")
lapply(p, require, character.only = TRUE)

enableWGCNAThreads()

#########################################
#  load data/inputs
#########################################

data_subset = "all" ## all nuclear wc ribo
home_dir = paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/")
# currently using counts based only on exons (for exons and introns see gene_from_intron_and_exon/ folder)
info = readRDS(paste0(home_dir, "info.RDS"))
vobj = readRDS(paste0(home_dir, "vobj.RDS"))
fit = readRDS(paste0(home_dir, "fit.RDS"))

# should be residuals
R = residuals(fit, vobj)
datExpr = as.matrix(t(R))
colnames(datExpr) = rownames(R)
rownames(datExpr) = colnames(R)

# prepare traits matrix
datTraits = model.matrix( ~ Cell_type + Method + gender, info) %>% as.data.frame

#########################################
#  Soft thresholding- choose correct beta
#########################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
library(WGCNA)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)



