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
p = c("annotate", "org.Mm.eg.db", "gplots", "variancePartition", "limma",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

data_subset = "all" ## all nuclear wc ribo
data_source = "all" ## exon all
home_dir = paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/")
# gene_from_intron_and_exon/
info = readRDS(paste0(home_dir, "info.RDS"))
# x = readRDS(paste0(home_dir, "/norm_from_", data_source, "_2018_07_23.RDS"))
# vobj_loc = paste0(home_dir, "/vobj_from_", data_source, "_2018_07_23.RDS")
varpart_loc = paste0(home_dir, "/varpart_from_", data_source, "_2018_08_28.RDS")

## confirm only male
info %>% group_by(gender) %>% tally

## create model matrix
# design = model.matrix(~ Cell_type + Method + 0, info) ## Cell_type + Method + gender

## estimate variance as function of mean expression
# vobj = voom(x, design, plot=TRUE)
# saveRDS(vobj, vobj_loc)
# vobj = readRDS(vobj_loc)
expr = readRDS("d1_d2_rnaseq/expression_data_fc/all/deseq2_from_all_vsd2018_08_28_ONLY_3WAY_GENES.RDS")
# expr = readRDS("d1_d2_rnaseq/expression_data_fc/all/deseq2_from_all_quantlog_2018_08_28.RDS")

form = ~ (1|Method) + (1|Cell_type) # + (1|gender)

## variance partition
# varPart = fitExtractVarPartModel(expr, form, info) ## vobj

# saveRDS(varPart, varpart_loc)
###
varPart = readRDS(varpart_loc)

## print to file (first sort by descending for Method)
sorted_vals = sort(varPart$Method, decreasing = T, index.return = T)$ix
out_df = sortCols(varPart[sorted_vals, ])
as.data.frame(out_df) %>% mutate(gene_ens_id = row.names(out_df)) %>%
  select(gene_ens_id, everything()) %>% 
  # filter(gene_ens_id %in% common_genes) %>% dim
  write_tsv("d1_d2_rnaseq/figures/variance_partition_2018_08_28/top_var_genes_2018_08_28.txt")
# top_variance_genes_2018_06_23.txt

# colnames(varPart) = c("Cell type", "Method", "Residuals") # "Gender",
# p = plotVarPart( sortCols(varPart[common_genes, ]), label.angle = 50)
# p
# ggsave(paste0("d1_d2_rnaseq/figures/variance_partition_fc_2018_06_18/var_explained.svg"),
#        p, width = 2.5, height = 3.25)

## looking at D1 and D2 expression across all datasets
drd1 = "ENSMUSG00000021478"
drd2 = "ENSMUSG00000032259"
xist = "ENSMUSG00000086503"
uty = "ENSMUSG00000068457"
lrrk2 = "ENSMUSG00000036273"


i = which.max(varPart$Method)
i = which(varPart$Cell_type > 0.7)
# mapIds(org.Mm.eg.db, row.names(varPart[i, ]), "SYMBOL","ENSEMBL")


#####################################
# Add gene names to VOBJ
#####################################

## for top genes, perform gene ontology enrichment
gene_name_vec = mapIds(org.Mm.eg.db, vobj$genes$GeneID, "SYMBOL","ENSEMBL")
vobj$genes$gene_name = ifelse(is.na(gene_name_vec), paste("Unk:",vobj$genes$GeneID), gene_name_vec)

row.names(info) = info$file_name
vobj$targets = bind_cols(vobj$targets, info)

#####################################
# heatmap of top genes for cell type
#####################################

topVarGenes = varPart %>% mutate(ens_name = row.names(varPart)) %>% arrange(desc(Cell_type)) %>% head(10)
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
color_scale = c(rep(color_scale[[1]], 40), color_scale, rep(color_scale[[length(color_scale)]], 40))
# 50 # 60

##
# len_kb = vobj$genes[ens_gene_list,"Length"]/1e3
# rpkm_norm_len = vobj$E[ens_gene_list, ]/len_kb

file_name = "d1_d2_rnaseq/figures/variance_partition_fc_2018_07_23/heatmap_cell_type.pdf"
# png(file_name, width = 6*300, height = 6*300, res = 300)
pdf(file_name, width = 6, height = 6)
p = heatmap.2(vobj$E[ens_gene_list,], # capped_expr, # vobj$E rpkm_norm_len, # 
          scale="row", ## none, row, column
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
          dendrogram="both")
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

i = lrrk2 #drd2 # drd2 drd1 lrrk2
len_kb = (vobj$genes[i,"Length"]/1e3)
## rpkm(x) cpm(x) vobj$E vobj$E/len_kb
GE = data.frame( Expression = vobj$E[i,], Cell_type = info$Cell_type, Method = info$Method)
# plotStratify( Expression ~ Method, GE) ## , colorBy=NULL, text=label, main=main
p = ggplot(GE, aes(x = Method, y = Expression, color = Cell_type)) + 
  geom_point() +
  scale_color_manual(values = c("red", "green")) + 
  theme_classic()
p
ggsave("d1_d2_rnaseq/figures/gene_expression_2018_07_23/lrrk2_expression_by_method_cell_type.png",
       p, width = 3, height = 3)

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
colList = c("royalblue", "red", "black") # for method
# order is: nuclear, ribo, wc
color.map = colList[num_label]

## label the gene expression level
color_scale = colorpanel(100, "Blue", "Black", "Yellow")[c(1:25, 45:55, 75:100)]
color_scale = c(rep(color_scale[[1]], 30), color_scale, rep(color_scale[[length(color_scale)]], 40))

file_name = "d1_d2_rnaseq/figures/variance_partition_fc_2018_07_23/heatmap_method.pdf"
# png(file_name, width = 6*300, height = 6*300, res = 300)
pdf(file_name, width = 6, height = 6)
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

