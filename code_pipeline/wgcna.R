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
p = c("WGCNA", "limma", #"coexpp",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr", "colorout")
lapply(p, require, character.only = TRUE)
source("d1_d2_rnaseq/code_pipeline/wgcna_plotting_functions.R")

enableWGCNAThreads(24)

#########################################
#  load data/inputs
#########################################

data_subset = "ribo" ## all nuclear wc ribo
home_dir = paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/")
# results_prefix = paste0(home_dir, "wgcna/vobjE_")
# vobj_scaled_noF_ vobjE_noF_

data_subset = "ribo" ## all nuclear wc ribo D1 D2 
# results_prefix = paste0("d1_d2_rnaseq/figures/wgcna_2018_04_18/", data_subset, "/vobjE_")
results_prefix = paste0("/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/", data_subset, "_results/vobjE_")

# currently using counts based only on exons (for exons and introns see gene_from_intron_and_exon/ folder)
info = readRDS(paste0(home_dir, "info.RDS"))
x = readRDS(paste0(home_dir, "norm.RDS"))

# re-define the design matrix, method, cell type, etc
# keep only a subset of samples
samples_to_keep = info %>% 
  # filter(gender == "M") %>% 
  filter(Method == data_subset) %>% ## disable for all
  # filter(Cell_type == data_subset) %>% 
  select(file_name) %>% unlist %>% as.character
length(samples_to_keep)
x = x[, samples_to_keep]
info %<>% filter(file_name %in% samples_to_keep)
# design = model.matrix(~ Cell_type + Method + gender + 0, info) #  + gender
# # design = model.matrix(~ Cell_type + 0, info) #  + gender
# vobj = voom(x, design, plot=TRUE) #
vobj = voom(x) # design, 

# alternatively use previously made voom object and fit that took cell type, method, and/or gender into account
# vobj = readRDS(paste0(home_dir, "vobj.RDS"))
# fit = readRDS(paste0(home_dir, "fit.RDS"))

# Only the residuals
# R = residuals(fit, vobj)
# datExpr = as.matrix(t(R))
# colnames(datExpr) = rownames(R)
# rownames(datExpr) = colnames(R)

# NOT residuals, NOT Gaussian
datExpr = as.matrix(t(vobj$E))
colnames(datExpr) = rownames(vobj$E)
rownames(datExpr) = colnames(vobj$E)

# NOT residuals, but normalized (subtract each row by rowMeans, then divide each row by row SD)
# vobj_mean = rowMeans(vobj$E)
# vobj_precision_sd = 1/sqrt(vobj$weights)
# vobj_centered = vobj$E %>% as.data.frame %>% mutate_all(funs(. - vobj_mean))
# vobj_scaled = as.matrix(vobj_centered)/as.matrix(vobj_precision_sd)
# datExpr = as.matrix(t(vobj_scaled))
# colnames(datExpr) = rownames(vobj_scaled)
# rownames(datExpr) = colnames(vobj_scaled)

# prepare traits matrix # + gender
datTraits = model.matrix( ~ Cell_type + Method + gender + 0, info) %>% as.data.frame %>% 
  mutate(MethodNuc = as.numeric((Methodribo == 0) & (Methodwc == 0))) %>% 
  select(Cell_typeD1:Methodwc, MethodNuc, genderM) #

datTraits = model.matrix( ~ Cell_type + gender + 0, info) %>% as.data.frame
# + gender + Method + Cell_type

#########################################
#  Soft thresholding- choose correct beta
#########################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2)) # , seq(from = 25, to=100, by=5)
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)
filename = paste0(results_prefix, "soft_thresholds_signed_18_04_19.pdf")
PlotSoftThreshold(sft, filename) 

#########################################
#  Coexpp library 
#########################################

beta_choice = 9
coexp_results_prefix = results_prefix %>% paste0(., "signed_beta", beta_choice, "_")
net = coexpressionAnalysis(datExpr, beta = beta_choice)

# 2415.987000 seconds
# ribo: 262.897000
# all: 3775.937000
# wc: 1400.099000


# Extract coexpp results into R from C++ object for saving
coexppClusters <- net$clusters
genePCTree <- net$genePCTree
geneModules <- net$geneModules
ims <- net$intraModularStatistics
rownames(ims) <- rownames(net$intraModularStatistics)
samples <- clusters(coexppClusters)$order

# plot dendrogram of genes with modules
pdf(file = paste0(coexp_results_prefix, "dendro_genes_coexpp_18_04_18.pdf"), width = 12, height = 9)
plotClustering(coexppClusters, geneModules, dendroLabels = FALSE)
dev.off()

# plot heatmap relationship between modules and traits
filename = paste0(coexp_results_prefix, "trait_module_18_04_18.pdf")
moduleTraitCor = plotTraitModule(datExpr, geneModules, datTraits, filename)

# plot correlation between traits
# filename = paste(results_prefix, "trait_cor_18_03_14.pdf", sep = "")
# plotTraitCor(datTraits, filename)

MEs = moduleEigengenes(datExpr, geneModules)$eigengenes

# createEigengeneNetwork = function(trait_interest, MEs, datTraits, coexp_results_prefix) {
# trait_interest = "genderM"
# # Cell_typeD2 Cell_typeD1 genderM MethodNuc Methodribo Methodwc
# trait.df = as.data.frame(datTraits[, trait_interest])
# names(trait.df) = trait_interest
trait_interest = "all"
trait.df = datTraits
MET = orderMEs(cbind(MEs, trait.df))
filename = paste0(coexp_results_prefix, "eigengene_dendro_", trait_interest, "_18_04_18.pdf")
pdf(file = filename, width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()
# }
# createEigengeneNetwork(trait_interest, MEs, datTraits, coexp_results_prefix)


# createGSMMTable(datExpr, geneModules, info.design.df, trait.interest, filename)


nSamples = nrow(datExpr)
MEsO = moduleEigengenes(datExpr, geneModules)$eigengenes
MEs = orderMEs(MEsO)
modNames = substring(names(MEs), 3)

trait_interest = "Cell_typeD1"
# Cell_typeD2 Cell_typeD1 genderM MethodNuc Methodribo Methodwc
# trait.df = as.data.frame(datTraits[, trait_interest])
# names(trait.df) = trait_interest

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(datTraits), sep="")
names(GSPvalue) = paste("p.GS.", names(datTraits), sep="")

# Create the starting data frame
geneInfo0 = data.frame(geneSymbol = colnames(datExpr),
                       moduleColor = geneModules,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for trait_interest
modOrder = order(-abs(cor(MEs, datTraits[, trait_interest], use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, 
                  -abs(geneInfo0[, paste("GS.", trait_interest, sep = "")]))
#geneOrder = order(geneInfo0$moduleColor, 
#                  -abs(geneInfo0[, paste("GS.", "L2HLHS.PHKG1", sep = "")]))
geneInfo = geneInfo0[geneOrder, ]

filename = paste0(coexp_results_prefix, "GS_MM_18_04_19.csv")
write.csv(geneInfo, file = filename)


# Export network to cytoscape format
# "midnightblue" yellow red
## for all, beta 16, unsigned:
modules = c("black", "yellow", "purple", "red") # 
# ribo beta 5:
modules = c("red", "blue", "brown", "yellow", "green")
probes = colnames(datExpr)
inModule = is.finite(match(geneModules, modules))
modProbes = probes[inModule]
# modTOM = 1 - tom.matrix[inModule, inModule]
# dimnames(modTOM) = list(modProbes, modProbes)
modTOM = dissTOM[modProbes, modProbes]

length(modTOM[modTOM > 0.01 & modTOM < 1])
length(modTOM[modTOM < 0.01])
#modTOM["BMP8B", modTOM["BMP8B", ] > 0.001]
modTOM["MTCH2", modTOM["MTCH2", ] > 0.01]
grep("MTCH2", rownames(modTOM))

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste0(coexp_results_prefix, "CytoscapeInput-edges-", 
                                                paste(modules, collapse="-"), ".5.txt"),
                               nodeFile = paste0(coexp_results_prefix, "CytoscapeInput-nodes-", 
                                                 paste(modules, collapse="-"), ".5.txt"),
                               weighted = TRUE,
                               # threshold = 0.01,
                               nodeNames = modProbes,
                               nodeAttr = geneModules[inModule])

## PLOT TOMs and adjacency matrices
## need to allocate 15 Gb of memory for these:
dissTOM <- tom(coexppClusters)[samples, samples]
# # note that this is 1-TOM, ie dissimilarity TOM
# geneTree = hclust(as.dist(dissTOM), method = "average");
# # Plot the resulting clustering tree (dendrogram)
# pdf(file = paste0(results_prefix, "beta10_dendro_genes_coexpp_TOM_18_03_14.pdf"), width = 12, height = 9)
# plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#      labels = FALSE, hang = 0.04)
# dev.off()
# adj.matrix <- adj(coexppClusters)[samples, samples]

######
#one-step

# coexp_file_base = "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/all_coexpp"
wgcna_file_base = paste0( "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/onestep_",
                          data_subset, "_signed_beta", beta_choice, "_")

net = blockwiseModules(datExpr, power = beta_choice,
                       TOMType = "signed", minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = wgcna_file_base,
                       verbose = 5)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
modMembers = data.frame(Gene = colnames(gnxp), Module = moduleColors)

# trying multi-step
adjacency = adjacency(gnxp, power = softPower)
str(adjacency)
TOM = TOMsimilarity(adjacency) #Topological Overlap Matrix
str(TOM)
dissTOM = 1-TOM
str(dissTOM)

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

