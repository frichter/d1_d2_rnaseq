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
p = c("WGCNA", "limma", "coexpp",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr", "colorout")
lapply(p, require, character.only = TRUE)
source("d1_d2_rnaseq/code_pipeline/wgcna_plotting_functions.R")

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
datTraits = model.matrix( ~ Cell_type + Method + gender + 0, info) %>% as.data.frame

#########################################
#  Soft thresholding- choose correct beta
#########################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
filename = paste0(home_dir, "wgcna/soft_thresholds_18_03_14.pdf")
PlotSoftThreshold(sft, filename) 

#########################################
#  Coexpp library 
#########################################

net = coexpressionAnalysis(datExpr, beta = 9)
# 2415.987000 seconds

coexp_file_base = "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/all_coexpp"

# Extract coexpp results into R from C++ object for saving
coexppClusters <- net$clusters
genePCTree <- net$genePCTree
geneModules <- net$geneModules
ims <- net$intraModularStatistics
rownames(ims) <- rownames(net$intraModularStatistics)
samples <- clusters(coexppClusters)$order
## need to allocate 15 Gb of memory for these:
# tom.matrix <- tom(coexppClusters)[samples, samples]
#note that this is actually 1-TOM
# adj.matrix <- adj(coexppClusters)[samples, samples]

# plot dendrogram of genes with modules
pdf(file = paste0(home_dir, "wgcna/dendro_genes_coexpp_18_03_14.pdf"), width = 12, height = 9)
plotClustering(coexppClusters, geneModules, dendroLabels = FALSE)
dev.off()

# plot heatmap relationship between modules and traits
filename = paste0(home_dir, "wgcna/trait_module_18_03_14.pdf")
moduleTraitCor = plotTraitModule(datExpr, geneModules, datTraits, filename)

# plot correlation between traits
filename = paste(home_dir, "wgcna/trait_cor_18_03_14.pdf", sep = "")
plotTraitCor(datTraits, filename)

# create MM and GS p-values table
# Mitral.Valve, 
gene_significance_tables = function(datExpr, geneModules, datTraits, home_dir, trait) {
  filename <- paste(home_dir, "wgcna/geneInfo_", trait, "_18_03_14.csv", sep = "")
  createGSMMTable(datExpr, geneModules, datTraits, trait, 
                  filename)
}

trait = names(datTraits)[[2]]
lapply(names(datTraits)[2:4], function(trait) 
  gene_significance_tables(datExpr, geneModules, datTraits, outdirectory, trait))

#########################################
#  Actually running WGCNA
#########################################


net = blockwiseModules(datExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = coexp_file_base, 
                       verbose = 3)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

map(net$blockGenes, length) %>% unlist %>% sum
map(net$dendrograms[[1]], head)
pdf(paste0(home_dir, "wgcna/geneDendro.pdf"), height=6, width=9)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms, mergedColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = paste0(coexp_file_base, "networkConstruction-auto.RData"))

#########################################
#  
#########################################

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(paste0(home_dir, "wgcna/module_trait_cor.pdf"), height=6, width=9)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()




