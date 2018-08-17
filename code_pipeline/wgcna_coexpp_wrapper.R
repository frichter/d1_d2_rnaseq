
#########################################
#  Coexpp library
#  WHY CAN I NOT SPECIFY SIGNED?? 
#########################################

beta_choice = 4
coexp_results_prefix = results_prefix %>% paste0(., "coexpp_beta", beta_choice, "_")
net = coexpressionAnalysis(datExpr, beta = beta_choice) ## how to pass the argument , TOMType = "signed"

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
pdf(file = paste0(coexp_results_prefix, "dendro_genes_coexpp_18_04_23.pdf"), width = 12, height = 9)
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
