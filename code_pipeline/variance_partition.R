# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 1/29/2018
# description: Variance partition
##############################################################

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("annotate", "org.Mm.eg.db", "gplots", "variancePartition",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

data_subset = "all" ## all nuclear wc ribo
home_dir = paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/")
# gene_from_intron_and_exon/
info = readRDS(paste0(home_dir, "info.RDS"))
vobj = readRDS(paste0(home_dir, "vobj.RDS"))

# remove female samples from analysis
files_to_keep = info %>% filter(gender == "M") %>% select(file_name) %>% unlist %>% as.character
vobj = vobj[, files_to_keep]
info %<>% filter(file_name %in% files_to_keep)


form = ~ (1|Method) + (1|Cell_type) # + (1|gender)

## variance partition
# varPart = fitExtractVarPartModel(vobj, form, info)
# saveRDS(varPart, "d1_d2_rnaseq/expression_data_fc/all/varpart_noF.RDS")
# varpart_noF.RDS  varpart.RDS
varPart = readRDS("d1_d2_rnaseq/expression_data_fc/all/varpart_noF.RDS")

## print to file (first sort by descending for Method)
# sorted_vals = sort(varPart$Method, decreasing = T, index.return = T)$ix
# out_df = sortCols(varPart[sorted_vals, ])
# as.data.frame(out_df) %>% mutate(gene_ens_id = row.names(out_df)) %>% 
#   select(gene_ens_id, everything()) %>% 
#   write_tsv("d1_d2_rnaseq/figures/variance_partition_fc_2018_06_18/top_variance_genes_2018_06_18.txt")

# colnames(varPart) = c("Cell type", "Method", "Residuals") # "Gender",
# p = plotVarPart( sortCols(varPart), label.angle = 50)
# p 
# ggsave(paste0("d1_d2_rnaseq/figures/variance_partition_fc_2018_06_18/var_explained.png"), 
#        p, width = 2.5, height = 3.25)

## looking at D1 and D2 expression across all datasets
drd1 = "ENSMUSG00000021478"
drd2 = "ENSMUSG00000032259"
xist = "ENSMUSG00000086503"
uty = "ENSMUSG00000068457"


i = which.max(varPart$Method)
i = which(varPart$Cell_type > 0.9)
## for top genes, perform gene ontology enrichment
gene_name_vec = mapIds(org.Mm.eg.db, vobj$genes$GeneID, "SYMBOL","ENSEMBL")
vobj$genes$gene_name = ifelse(is.na(gene_name_vec), paste("Unk:",vobj$genes$GeneID), gene_name_vec)

row.names(info) = info$file_name
vobj$targets = bind_cols(vobj$targets, info)

################################
# heatmap of top genes for gender
################################
topVarGenes = varPart %>% mutate(ens_name = row.names(varPart)) %>% arrange(desc(gender)) %>% head(7)
topVarGenes %<>% 
  mutate(gene_name = mapIds(org.Mm.eg.db, topVarGenes$ens_name, "SYMBOL","ENSEMBL")) %>% 
  ## for top 7 gender genes use:
  mutate(gene_name = ifelse(is.na(gene_name), "Gm26992", gene_name))
topVarGenes

ens_gene_list = topVarGenes$ens_name
# mycol <- colorpanel(1000,"blue","white","red")
file_name = "d1_d2_rnaseq/figures/variance_partition_fc_2018_01_29/heatmap_gender.png"
png(file_name, width = 6*300, height = 5*300, res = 300)
heatmap.2(vobj$E[ens_gene_list,], 
          scale="row",
          labRow = vobj$genes[ens_gene_list,"gene_name"],
          labCol="", ## vobj$targets$gender
          trace="none", density.info="none", 
          margin=c(12,18),
          # lhei=c(2,5),
          # lwid = c(1,4),
          # scale="none",
          # ColSideColors = color.map, 
          col = colorpanel(50, "Blue", "Black", "Yellow"), 
          # cexRow = 1, 
          # cexCol = 0.4, 
          key = FALSE, symkey=FALSE,
          dendrogram="column")
dev.off()


#####################################
# heatmap of top genes for cell type
#####################################

topVarGenes = varPart %>% mutate(ens_name = row.names(varPart)) %>% arrange(desc(Cell_type)) %>% head(20)
topVarGenes %<>% 
  mutate(gene_name = mapIds(org.Mm.eg.db, topVarGenes$ens_name, "SYMBOL","ENSEMBL")) 
topVarGenes

ens_gene_list = topVarGenes$ens_name

## label the x-axis covariate (e.g., Cell type, method)
num_label = as.numeric(vobj$targets$Cell_type) # Method Cell_type
# colList = rainbow(unique(num_label) + 1)
colList = c("Green", "Red") # for cell type
# colList = c("purple3", "skyblue", "red3") # for method
color.map = colList[num_label]

## label the gene expression level
color_scale = colorpanel(100, "Blue", "Black", "Yellow")[c(1:25, 45:55, 75:100)]
color_scale = c(rep(color_scale[[1]], 50), color_scale, rep(color_scale[[length(color_scale)]], 60))

file_name = "d1_d2_rnaseq/figures/variance_partition_fc_2018_06_18/heatmap_cell_type.png"
png(file_name, width = 6*300, height = 6*300, res = 300)
p = heatmap.2(vobj$E[ens_gene_list,], # capped_expr, # 
          scale="row",
          labRow = vobj$genes[ens_gene_list,"gene_name"],
          labCol="", # vobj$targets$Cell_type,
          trace="none", density.info="none", 
          margin=c(5,16),
          # lhei=c(2,5),
          # lwid = c(1,4),
          # scale="none",
          ColSideColors = color.map,
          col = color_scale, # colorpanel(100, "Blue", "Black", "Yellow")[c(1:25, 45:55, 75:100)],
          # cexRow = 1, 
          # cexCol = 0.4, 
          # key = FALSE, 
          # symkey=FALSE,
          dendrogram="row")
dev.off()

drd1_sig = topVarGenes[rev(p$rowInd)[c(4:9)], ]
drd1_sig$gene_name %>% paste(collapse = ", ")
## D1 cluster: Sfxn1, Cntnap3, Ttc23, Drd1, Pdyn, Dlk1

drd2_sig = topVarGenes[rev(p$rowInd)[c(c(1:3),c(10:20))], ]
drd2_sig$gene_name %>% paste(collapse = ", ")
## cluster 1: Gpr6, P2ry1, Sp9
## cluster 2: Drd2, Adora2a, NA, Gnb5, Oprd1, NA, Ankk1, Ndnf, Ttc12, Penk, Gpr88
## consistent with this, cluster 2 most enriched GO terms are 
# chemical synaptic transmission, postsynaptic (GO:0099565) and synaptic transmission, dopaminergic (GO:0001963) (both 8.334e-05) 

# # i = which(varPart$Method > 0.90)
# rownames(varPart)[i]
# i = c(xist, uty) ## c(drd1, drd2)
# plotVarPart( sortCols(varPart[t30$gene_name, ]), label.angle=50 , main=rownames(varPart[i,])) ##
# 
# data.frame(t(vobj$E[t30$gene_name, ]), Cell_type = info$Cell_type) %>% group_by(Cell_type) %>% 
#   summarise_all(mean) %>% ungroup %>% as.data.frame
# 
# ensIDs = row.names(vobj$E[i, ]) 
# 
# GE = data.frame( Expression = vobj$E[drd2,], Cell_type = info$Cell_type, Method = info$Method)
# plotStratify( Expression ~ Cell_type, GE) ## , colorBy=NULL, text=label, main=main


#####################################
# heatmap of top genes for method
#####################################

topVarGenes = varPart %>% mutate(ens_name = row.names(varPart)) %>% arrange(desc(Method)) %>% head(40)
topVarGenes %<>% 
  mutate(gene_name = mapIds(org.Mm.eg.db, topVarGenes$ens_name, "SYMBOL","ENSEMBL")) 
topVarGenes

ens_gene_list = topVarGenes$ens_name

## label the x-axis covariate (e.g., Cell type, method)
num_label = as.numeric(vobj$targets$Method) # Method Cell_type
# colList = c("purple3", "skyblue", "red3") # for method
colList = c("blue", "red", "black") # for method
# order is: nuclear, ribo, wc
color.map = colList[num_label]

## label the gene expression level
color_scale = colorpanel(100, "Blue", "Black", "Yellow")[c(1:25, 45:55, 75:100)]
color_scale = c(rep(color_scale[[1]], 30), color_scale, rep(color_scale[[length(color_scale)]], 40))

file_name = "d1_d2_rnaseq/figures/variance_partition_fc_2018_06_18/heatmap_method.png"
png(file_name, width = 6*300, height = 6*300, res = 300)
p = heatmap.2(vobj$E[ens_gene_list,], 
              scale="row",
              labRow = vobj$genes[ens_gene_list,"gene_name"],
              labCol="", # vobj$targets$Method,
              trace="none", density.info="none", 
              margin=c(5,16),
              # lhei=c(2,5),
              # lwid = c(1,4),
              # scale="none",
              ColSideColors = color.map,
              col = color_scale,
              # cexRow = 1, 
              # cexCol = 0.4, 
              # key = FALSE, 
              # symkey=FALSE,
              dendrogram="row")
dev.off()

