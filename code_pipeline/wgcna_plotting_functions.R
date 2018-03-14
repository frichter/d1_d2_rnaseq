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


plotTraitModule = function(datExpr, geneModules, datTraits, filename) {
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
  
  pdf(file = filename, width = 10, height = 6)
  par(mar = c(6, 8.5, 3, 3))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(datTraits),
                 yLabels = colnames(MEs),
                 ySymbols = colnames(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.2,
                 cex.lab.y = 0.5, 
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
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
                 xLabels = colnames(info.design.df),
                 yLabels = colnames(info.design.df),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.3,
                 cex.lab.y = 0.5, 
                 cex.lab.x = 0.5, 
                 zlim = c(-1,1),
                 main = paste("Trait-Trait Relationships"))
  dev.off()
}


