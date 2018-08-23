# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 3/13/2018
# description: helper functions for wcgna.R
##############################################################


PlotSoftThreshold = function(sft, filename) {
  pdf(filename, height=6, width=14)
  par(mfrow = c(1,2))
  cex1 = 0.9
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.85,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}


plotTraitModule = function(datExpr, geneModules, datTraits, filename, color_scale,
                           width = 4, height = 6) {
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  MEsO = moduleEigengenes(datExpr, geneModules)$eigengenes
  MEs = orderMEs(MEsO)
  
  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  
  pdf(file = filename, width = width, height = height)
  par(mar = c(6, 8.5, 3, 3))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(datTraits),
                 yLabels = colnames(MEs),
                 ySymbols = colnames(MEs),
                 colorLabels = FALSE,
                 colors = color_scale,
                 # colors = blueWhiteRed(50),
                 # cexCol = "#FFFFFF", ## not a real argument, unclear how to specify the textMatrix colors
                 # textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.2,
                 cex.lab.y = 0.5, 
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
  mtx_filename = gsub(".pdf", ".txt", filename)
  
  mod_trait_cor = cbind(signif(moduleTraitCor, 2), signif(moduleTraitPvalue, 1)) %>% as.data.frame
  traits = ncol(datTraits)
  names(mod_trait_cor) = paste0(names(mod_trait_cor), c(rep("cor", traits), rep("p", traits)))
  mod_trait_cor %>% 
    mutate(module = row.names(mod_trait_cor)) %>% 
    select(module, everything()) %>% 
    write_tsv(mtx_filename)
  return(moduleTraitCor)
}

plotTraitCor = function(datTraits, filename) {
  traitCor = cor(datTraits, datTraits, use = "p")
  
  # Will display correlations and their p-values
  textMatrix =  paste(signif(traitCor, 2), sep = "")
  dim(textMatrix) = dim(traitCor)
  
  pdf(file = filename, width = 10, height = 6)
  par(mar = c(6, 8.5, 3, 3))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = traitCor,
                 xLabels = colnames(datTraits),
                 yLabels = colnames(datTraits),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.3,
                 cex.lab.y = 0.5, 
                 cex.lab.x = 0.5, 
                 zlim = c(-1,1),
                 main = paste("Trait-Trait Relationships"))
  dev.off()
}


createGSMMTable = function(datExpr, moduleColors, datTraits, 
                            trait_interest, filename) {
  
  nSamples = nrow(datExpr)
  MEsO = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEsO)
  modNames = substring(names(MEs), 3)
  
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
                         moduleColor = moduleColors,
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
                    -abs(geneInfo0[, paste("GS.", make.names(trait_interest), sep = "")]))
  #geneOrder = order(geneInfo0$moduleColor, 
  #                  -abs(geneInfo0[, paste("GS.", "L2HLHS.PHKG1", sep = "")]))
  geneInfo = geneInfo0[geneOrder, ]
  
  write.csv(geneInfo, file = filename)
}



# createModuleScatterPlots <- function(datExpr, info.design.df, geneModules, 
#                                      geneModules.interest, filename) {
#   MEsO = moduleEigengenes(datExpr, geneModules)$eigengenes
#   MEs = orderMEs(MEsO)
#   modNames = substring(names(MEs), 3)
#   
#   trait.df = as.data.frame(info.design.df[, trait.interest])
#   names(trait.df) = trait.interest
#   geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
#   geneTraitSignificance = as.data.frame(cor(datExpr, trait.df, use = "p"))
#   names(geneTraitSignificance) = paste("GS.", names(trait.df), sep="")
#   
#   pdf(file = filename, width = 7, height = 7)
#   for(module in geneModules.interest) {
#     column = match(module, modNames)
#     moduleGenes = geneModules==module
#     par(mfrow = c(1,1))
#     verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                        abs(geneTraitSignificance[moduleGenes, 1]),
#                        xlab = paste("Module Membership in", module, "module"),
#                        ylab = paste("Gene significance for", trait.interest), 
#                        main = paste("Module membership vs. gene significance\n"),
#                        cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#   }
#   dev.off()
# }


WrapperForCytoscapeExport = function(modules, datExpr, TOM, moduleColors, cytoprefix) {
  print(modules)
  # Select module probes
  genes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules))
  modGenes = genes[inModule]
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modGenes, modGenes)
  print(modTOM[1:5, 1:5])
  print(dim(modTOM))
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste0(cytoprefix, "/cs_edges_",
                                                   paste(modules, collapse="-"), ".txt"),
                                 nodeFile = paste0(cytoprefix, "/cs_nodes_",
                                                   paste(modules, collapse="-"), ".txt"),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
  return(modules)
}

# createEigengeneNetwork <- function(datExpr, geneModules, trait.interest, filename){
#   MEs = moduleEigengenes(datExpr, geneModules)$eigengenes
#   trait.df = as.data.frame(info.design.df[, trait.interest])
#   names(trait.df) = trait.interest
#   MET = orderMEs(cbind(MEs, trait.df))
#   pdf(file = filename, width = 6, height = 6)
#   par(cex = 1.0)
#   plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
#                         plotHeatmaps = FALSE)
#   par(cex = 1.0)
#   plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
#                         plotDendrograms = FALSE, xLabelsAngle = 90)
#   dev.off()
# }
# 
# 
# HeatmapWrapper <- function(m, colors, plot.raise_power=4, ...) {
#   diag(m) = NA
#   m = m ^ plot.raise_power
#   
#   heatmap(
#     m,
#     Rowv=NA, Colv=NA, scale="none", revC=TRUE, symm=TRUE, labRow="", labCol="",     
#     ColSideColors=as.character(colors),
#     RowSideColors=as.character(colors),
#     ...
#   )
# }
