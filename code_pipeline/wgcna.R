# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 3/13/2018
# description: WGCNA
##############################################################

# export TMP=/sc/orga/projects/chdiTrios/Felix/dna_rna/eqtl_wgs/fastqtl_2018_01/
# module load R/3.3.1 ## use this minerva version for for coexpp ## R/3.4.1 ##
# R

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
setwd("/hpc/users/richtf01/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("WGCNA", "limma", "gplots", #"coexpp",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr", "colorout")
lapply(p, require, character.only = TRUE)
source("d1_d2_rnaseq/code_pipeline/wgcna_plotting_functions.R")

enableWGCNAThreads(6)

#########################################
#  load data/inputs
#########################################

parent_subset = "ribo" ## all nuclear wc ribo
home_dir = paste0("d1_d2_rnaseq/expression_data_fc/", parent_subset, "/")
# results_prefix = paste0(home_dir, "wgcna/vobjE_")
# vobj_scaled_noF_ vobjE_noF_

data_subset = "ribo" ## all nuclear wc ribo D1 D2 
# results_prefix = paste0("/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/", data_subset, "_results/vobjE_")
results_prefix = paste0("/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/", data_subset, "_results/vst_")

info = readRDS(paste0(home_dir, "info.RDS"))

# keep only a subset of samples (only need to run this for D1 and D2)
# samples_to_keep = info %>% 
#   # filter(gender == "M") %>% 
#   # filter(Method == data_subset) %>% ## disable for all
#   filter(Cell_type == data_subset) %>%
#   select(file_name) %>% unlist %>% as.character
# length(samples_to_keep)
# info %<>% filter(file_name %in% samples_to_keep)

# prepare traits matrix # + gender
datTraits = model.matrix( ~ Cell_type + Method + gender + 0, info) %>% as.data.frame %>% 
  mutate(MethodNuc = as.numeric((Methodribo == 0) & (Methodwc == 0))) %>% 
  select(Cell_typeD1:Methodwc, MethodNuc, genderM) #

datTraits = model.matrix( ~ Cell_type + gender + 0, info) %>% as.data.frame
# + gender + Method + Cell_type
head(datTraits)
names(datTraits) = c("D1 neurons ", "D2 neurons ", "Sex (Male) ")

#############################################
#  Using the voom-limma log-transformed data:
#############################################

# currently using counts based only on exons (for exons and introns see gene_from_intron_and_exon/ folder)
x = readRDS(paste0(home_dir, "norm.RDS"))
x = x[, samples_to_keep]

# design = model.matrix(~ Cell_type + Method + gender + 0, info) #  + gender
# vobj = voom(x, design, plot=TRUE) #
vobj = voom(x) # design, 

# NOT residuals, NOT Gaussian
datExpr = as.matrix(t(vobj$E))
colnames(datExpr) = rownames(vobj$E)
rownames(datExpr) = colnames(vobj$E)

#############################################
#  Alternatively variance stabilized DESeq2
#############################################

# Using the DESeq2 input
vsd = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", parent_subset, "/deseq2_vsd2018_04_23.RDS"))
## only subset for D1 and D2 (other ones are already subset)
# vsd = vsd[, samples_to_keep]

datExpr = as.matrix(t(vsd))
colnames(datExpr) = rownames(vsd)
rownames(datExpr) = colnames(vsd)

#########################################
#  Soft thresholding- choose correct beta
#########################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2)) # , seq(from = 25, to=100, by=5)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)
filename = paste0(results_prefix, "soft_thresholds_signed_18_04_26.pdf")
PlotSoftThreshold(sft, filename) 

################################################
# one-step pure WGCNA (ie no coexpp library)
################################################

beta_choice = 5
wgcna_file_base = results_prefix %>% paste0(., "bicor_signed_beta", beta_choice, 
                                            "_min100_mergecutheight2neg2_static99_",
                                            "minKMEtoStay1neg2_pamF_")
# max power for signed is 30..
net = blockwiseModules(datExpr, power = beta_choice,
                       networkType = "signed",
                       TOMType = "signed", 
                       detectCutHeight = 0.99,
                       minModuleSize = 100,
                       reassignThreshold = 0,
                       minKMEtoStay = 0.1, 
                       mergeCutHeight = 0.2, # (1 - correlation btw eigengenes) for merging modules
                       corType="bicor",
                       numericLabels = TRUE,
                       pamStage = FALSE, 
                       # this also greatly determines how "clean". Geschwind uses negative (probably means false)
                       # pamRespectsDendro = TRUE,
                       maxBlockSize = 46300,
                       saveTOMs = TRUE,
                       saveTOMFileBase = wgcna_file_base,
                       verbose = 3)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
modMembers = data.frame(Gene = colnames(datExpr), Module = moduleColors)

unique(moduleColors) %>% length

saveRDS(net, paste0(wgcna_file_base, "bwm_out_18_05_05.RDS"))

# Plot the dendrogram and the module colors underneath
# pdf(file = paste0(wgcna_file_base, "dendro_genes_min20_18_04_23.pdf"), width = 12, height = 9)
pdf(file = paste0(wgcna_file_base, "dendro_genes_18_05_05.pdf"),
    width = 9, height = 4.5)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# get module eigengenes
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

## label the gene expression level
# library(gplots)
color_scale = colorpanel(50, "Blue", "Black", "Yellow")
# displayColors(color_scale)

MET = orderMEs(cbind(MEs, datTraits))
filename = paste0(wgcna_file_base, "eigengene_dendro_18_05_05.pdf")
pdf(file = filename, width = 6, height = 6)
# Plot the relationships among the eigengenes and the trait 
# Plot the dendrogram
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      signed = TRUE,
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap",
                      # marHeatmap = c(3,4,2,2),# previous
                      marHeatmap = c(7,7,2,2),
                      signed = TRUE,
                      heatmapColors = color_scale,
                      plotDendrograms = FALSE,
                      xLabelsAngle = 90)
dev.off()

# plot heatmap relationship between modules and traits
filename = paste0(wgcna_file_base, "trait_module_18_05_05.pdf")
moduleTraitCor = plotTraitModule(datExpr, moduleColors, datTraits, filename, color_scale)

# sending gene modules to files

trait_interest = "D1 neurons "
# Cell_typeD2 Cell_typeD1 genderM Methodnuclear Methodribo Methodwc
filename = paste0(wgcna_file_base, "GS_MM_18_05_05.csv")
createGSMMTable(datExpr, moduleColors, datTraits, trait_interest, filename)


####################################
# Export every module to cytoscape
####################################
load(paste0(wgcna_file_base, "-block.1.RData"))

TOM = TOMsimilarityFromExpr(datExpr, power = beta_choice)

# Select modules
modules = unique(moduleColors)[[2]]
cytoprefix = gsub("vst_", "", results_prefix)
WrapperForCytoscapeExport(modules, datExpr, TOM, moduleColors, cytoprefix) 


## link for GO enrichment:
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-04-Interfacing.R



# trying multi-step
# adjacency = adjacency(gnxp, power = softPower)
# str(adjacency)
# TOM = TOMsimilarity(adjacency) #Topological Overlap Matrix
# str(TOM)
# dissTOM = 1-TOM
# str(dissTOM)

####################################
# Previous coexpp code to parse:
####################################

# create MM and GS p-values table
# Mitral.Valve, 
# gene_significance_tables = function(datExpr, geneModules, datTraits, home_dir, trait) {
#   filename <- paste(home_dir, "wgcna/vobjE_geneInfo_", trait, "_18_03_14.csv", sep = "")
#   createGSMMTable(datExpr, geneModules, datTraits, trait, 
#                   filename)
# }
# 
# trait = names(datTraits)[[2]]
# lapply(names(datTraits)[2:4], function(trait) 
#   gene_significance_tables(datExpr, geneModules, datTraits, outdirectory, trait))
# 
# # create MM and GS p-values table
# filename <- paste(outdirectory, "geneInfo.csv", sep = "")
# createGSMMTable(datExpr, geneModules, info.design.df, trait.interest, 
#                 filename)
# 
# # create scatter plots of GS and MM for every gene in a module
# filename <- paste(outdirectory, "module_scatter_plots_", trait.interest, ".pdf", sep = "")
# geneModules.interest <- c("red", "black", "midnightblue", "peru", "gold", 
#                           "sienna", "magenta", "lightgreen", "salmon", "yellow", "coral")
# #geneModules.interest <- c("tan", "gold", "sienna")
# #geneModules.interest <- c("peru", "sienna", "magenta", "lightgreen", 
# #    "salmon", "coral")
# createModuleScatterPlots(datExpr, info.design.df, geneModules, 
#                          geneModules.interest, filename)
# 
# # create eigengene dendrogram and network to identify closely related modules
# filename <- paste(outdirectory, "eigengene_network_", trait.interest, ".pdf", sep = "")
# createEigengeneNetwork(datExpr, geneModules, trait.interest, filename)
# 
# # create a topological overlap matrix plot. Must use original coexpp C++ object! cannot use R tom.matrix
# filename <- paste(outdirectory, "plot_TOM.pdf", sep = "")
# pdf(filename)
# samplingThreshold=1000
# samples <- clusters(coexppClusters)$order[sort(sample(length(coexppClusters), samplingThreshold))]
# HeatmapWrapper(
#   sampleMatrix(coexppClusters, samples, kind="tom"),
#   geneModules[samples], plot.raise_power = 10
# )
# dev.off()
# 

