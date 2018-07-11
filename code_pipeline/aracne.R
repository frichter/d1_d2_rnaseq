
# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 6/10/2018
# description: Run aracne
##############################################################

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("limma", "edgeR", "annotate", "org.Mm.eg.db", "DESeq2", "minet", "topGO", "goseq",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)


data_subset = "nuclear" ## all nuclear wc ribo
## load module results
wgcna_dir = paste0("d1_d2_rnaseq/figures/wgcna_2018_04_18/", data_subset, "/cytoscape_modules")
module_list = list.files(wgcna_dir, "cs_edges.*txt", full.names = T)
# which files correspond the blue/salmon modules
module_color = "turquoise" # blue brown grey turquoise
grep(module_color, module_list)

### running ARACNE on the WGCNA edges
wgcna_network = read_tsv(module_list[[15]])


from_nodes = wgcna_network$fromNode %>% unique
to_nodes = wgcna_network$toNode %>% unique
nodes_to_keep = c(from_nodes, to_nodes) %>% unique
# note that there are typically nodes unique to fromNode and toNode (typically 1 but sometimes more)
length(to_nodes)
length(nodes_to_keep)


### load same expression data used for WGCNA
vsd = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/deseq2_vsd2018_04_23.RDS"))

# subset the expression to just genes of interest
# confirm all nodes being kept have expression data (bumbers should be the same)
length(nodes_to_keep)
sum(row.names(vsd) %in% nodes_to_keep)
module_expr = vsd[row.names(vsd) %in% nodes_to_keep, ] %>% t %>% as.data.frame
# calculate the mutual information criteria. Note that MI methods need discretized data
# calculate rank
# module_expr_rank = matrix(rank(module_expr, ties = "min"), ncol = ncol(module_expr))
# same results for ranked and original data
module_mim = build.mim(module_expr, estimator = "mi.empirical",
                       disc = "equalfreq")

# library(infotheo)
# ?mutinformation
# mutinformation

module_aracne = aracne(module_mim)

# compare the number of edges with ARACNE and WGCNA
sum(module_aracne != 0)
dim(wgcna_network)

# convert to long format for cytoscape import
module_aracne_long = module_aracne %>% as.data.frame %>% 
  mutate(fromNode = row.names(module_aracne)) %>% 
  gather(key = "toNode", value = "mic", -fromNode) %>% 
  # remove edges with mic==0
  filter(mic != 0)

hist(module_aracne_long$mic)

# join with WGCNA data
module_aracne_long %<>% 
  arrange(fromNode) %>% 
  mutate(direction = "undirected") %>% 
  left_join(wgcna_network %>% select(fromNode:weight)) %>% 
  rename(mic_aracne = mic, weight_wgcna = weight) 

# remove edges with NA according to WGCNA
module_aracne_wgcna = module_aracne_long %>% filter(!is.na(weight_wgcna))
dim(module_aracne_wgcna)

# write to file
module_aracne_wgcna %>% 
  write_tsv(paste0(wgcna_dir, "/aracne_fromR_", module_color, ".txt"))

module_aracne_wgcna %>% 
  select(-weight_wgcna, -mic_aracne, -direction) %>% 
  write_tsv(paste0(wgcna_dir, "/aracne_fromR_", module_color, "_noWeights.txt"))

###########################
# load aracne/wgcna module
# to decrease size
###########################

module_aracne_wgcna = read_tsv(paste0(wgcna_dir, "/aracne_fromR_", module_color, ".txt"))

# any other filters to decrease module size?
module_aracne_wgcna %>% filter(!is.na(weight_wgcna)) %>% dim
# ggplot(module_aracne_wgcna, aes(x = mic_aracne)) + geom_density() + theme_classic()
# ggplot(module_aracne_wgcna, aes(x = weight_wgcna)) + geom_density() + theme_classic()
# ggplot(module_aracne_wgcna, aes(x = weight_wgcna, y = mic_aracne)) + geom_point() + theme_classic()

# remove the bottom quartile of edges
module_aracne_wgcna_top75 = module_aracne_wgcna %>% 
  filter(weight_wgcna > quantile(weight_wgcna, 0.25),
         mic_aracne > quantile(mic_aracne, 0.25))

module_aracne_wgcna %>% #head(100000) %>% 
  select(weight_wgcna, mic_aracne) %>% cor

# write to file (without weights to minimize size)
module_aracne_wgcna_top75 %>% 
  select(-weight_wgcna, -mic_aracne, -direction) %>% 
  write_tsv(paste0(wgcna_dir, "/aracne_fromR_", module_color, "_t75_noWeights.txt"))

##################################
## GO enrichment
##################################

# create dataframe containing all genes with those in ARACNE module having a nonzero_mic=
# This way only testing enrichment among genes where expression was profiled
final_aracne_nodes = c(module_aracne_wgcna$fromNode, module_aracne_wgcna$toNode) %>% unique

expr_genes = unique(row.names(vsd))
all_genes_w_scores = as.numeric(expr_genes %in% final_aracne_nodes)
sum(all_genes_w_scores)
length(all_genes_w_scores)
names(all_genes_w_scores) = expr_genes

# get GO IDs from gene names
gene_map = getgo(names(all_genes_w_scores), 'mm9', 'ensGene')
# look for NAs in the gene_map (how many genes are represented on GO?)
sum(is.na(names(gene_map)))
sum(!is.na(names(gene_map)))
gene_map = gene_map[!is.na(names(gene_map))]

# function to return indices of genes in ARACNE module
mySelGenes = function(score) {
  return (score != 0)
}

GetTopTermsPerOntology = function(ontology_i, all_genes_w_scores, gene_map) {
  print(ontology_i)
  module_go = new("topGOdata",
                  description = "GO terms associated with an ARACNE module",
                  ontology = ontology_i, # CC MF BP
                  allGenes = all_genes_w_scores, 
                  geneSel = mySelGenes,
                  annot = annFUN.gene2GO,
                  gene2GO = gene_map)
  
  test_statistic = new("classicCount", testStatistic = GOFisherTest, name = "FET")
  results_fet = getSigGroups(module_go, test_statistic)
  print(results_fet)
  results_final = GenTable(module_go, classic = results_fet, topNodes = 10)
  results_final %<>% mutate(ontology = ontology_i)
  return(results_final)
}

go_results = map_df(c("CC", "MF", "BP"), GetTopTermsPerOntology, all_genes_w_scores, gene_map)
go_results
write_tsv(go_results, paste0(wgcna_dir, "/aracne_", module_color, "_GO.txt"))

