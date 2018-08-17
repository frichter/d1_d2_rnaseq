
# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 6/10/2018
# description: Run aracne
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
p = c("limma", "edgeR", "annotate", "org.Mm.eg.db", "DESeq2", "minet", "topGO", "goseq", "colorout",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)


data_subset = "ribo_w_Female" ## all nuclear wc ribo ribo_w_Female
## load module results
wgcna_dir = paste0("/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/", 
                   data_subset, "_results/cytoscape_modules")
wgcna_dir = paste0("d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/",
                   data_subset, "/cytoscape_modules")
module_list = list.files(wgcna_dir, "cs_edges.*txt", full.names = T)

### load same expression data used for WGCNA
vsd = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset,
                     "/deseq2_from_all_vsd2018_07_12.RDS"))

# which files correspond the blue/salmon modules
module_color = "grey" # blue brown grey turquoise
grep(paste0("_", module_color, ".txt"), module_list)

### running ARACNE on the WGCNA edges
module_i = module_list[[5]]

Generate_aracne_module = function(module_i, vsd, wgcna_dir) {
  module_color = gsub(".*_edges_", "", module_i) %>% gsub(".txt", "", .)
  print(module_color)
  
  aracne_out_file = paste0(wgcna_dir, "/aracne_", module_color, ".txt")
  if(file.exists(aracne_out_file)) {
    print(paste(module_color, "already ran through aracne"))
    return(module_color)
  }
  

  wgcna_network = read_tsv(module_i)
  ## for large modules, read nodes instead of WGCNA network!!
  # wgcna_nodes = read_tsv(gsub("_edges", "_nodes", module_i))
  
  # remove the bottom quantile of edges
  # wgcna_network %<>% filter(weight > quantile(weight, 0.5))
  
  ## get list of genes
  from_nodes = wgcna_network$fromNode %>% unique
  to_nodes = wgcna_network$toNode %>% unique
  
  nodes_to_keep = c(from_nodes, to_nodes) %>% unique
  print("Nodes being kept, confirm they are all about the same length:")
  # note that there are typically nodes unique to fromNode and toNode (typically 1 but sometimes more)
  print(length(to_nodes))
  print(length(nodes_to_keep))
  
  # subset the expression to just genes of interest
  # confirm all nodes being kept have expression data (numbers should be the same)
  # print(length(nodes_to_keep))
  print(sum(row.names(vsd) %in% nodes_to_keep))
  module_expr = vsd[row.names(vsd) %in% nodes_to_keep, ] %>% t %>% as.data.frame
  ## FOR LARGE MODULES:
  # module_expr = vsd[row.names(vsd) %in% wgcna_nodes$nodeName, ] %>% t %>% as.data.frame
  
  # calculate the mutual information criteria. Note that MI methods need discretized data
  # calculate rank
  # module_expr_rank = matrix(rank(module_expr, ties = "min"), ncol = ncol(module_expr))
  # same results for ranked and original data
  module_mim = build.mim(module_expr, estimator = "mi.empirical",
                         disc = "equalfreq")
  module_aracne = aracne(module_mim)
  
  # saveRDS(module_aracne, gsub(".txt", "_temp_aracne_wide.RDS", aracne_out_file))

  # compare the number of edges with ARACNE and WGCNA
  print("Number of non-zero edges in aracne and WGCNA")
  print(sum(module_aracne != 0))
  print(dim(wgcna_network))
  
  print("convert to long format for cytoscape import..")
  module_aracne_long = module_aracne %>% as.data.frame %>% 
    mutate(fromNode = row.names(module_aracne)) %>% 
    gather(key = "toNode", value = "mic", -fromNode) %>% 
    # remove edges with mic==0
    filter(mic != 0)
  
  # module_aracne_long %>%
  #   write_tsv(gsub("aracne_", "temp_intermediate_aracne_long_", aracne_out_file))

  # hist(module_aracne_long$mic)
  # join with WGCNA data
  module_aracne_long %<>% 
    arrange(fromNode) %>% 
    mutate(direction = "undirected") %>% 
    left_join(wgcna_network %>% select(fromNode:weight)) %>% 
    rename(mic_aracne = mic, weight_wgcna = weight) 
  # remove edges with NA according to WGCNA
  module_aracne_wgcna = module_aracne_long %>% filter(!is.na(weight_wgcna))
  print(dim(module_aracne_wgcna))
  # write to file
  module_aracne_wgcna %>% 
    write_tsv(aracne_out_file)
  return(module_color)
}

map(module_list, Generate_aracne_module, vsd, wgcna_dir)


#############################
# Get all edges for top 20
# nodes (from kME)
#############################

## local directories:
gs_mm_loc_all = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/all/vst_bicor_signed_beta16_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"
gs_mm_loc_nuc = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/nuclear/vst_bicor_signed_beta18_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"
gs_mm_loc_wc = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/wc/vst_bicor_signed_beta14_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"
gs_mm_loc_ribo = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/ribo/vst_bicor_signed_beta9_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"
gs_mm_loc_ribo_w_F = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/ribo_w_Female/vst_bicor_signed_beta9_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"

gs_mm_loc_list = list(gs_mm_loc_all, gs_mm_loc_nuc, gs_mm_loc_wc, gs_mm_loc_ribo, gs_mm_loc_ribo_w_F)
names(gs_mm_loc_list) = c("all", "nuclear", "wc", "ribo" , "ribo_w_Female")

data_subset = "ribo_w_Female"
gene_mm_gs = read_csv(gs_mm_loc_list[[data_subset]])

gene_mm_gs %<>% select(-X1)

## add _ to end so you can match to end of line
names(gene_mm_gs)[7:ncol(gene_mm_gs)] = paste0(names(gene_mm_gs)[7:ncol(gene_mm_gs)], "_")

# module_color_i = module_list[[1]]

get_edges_for_t20_nodes = function(module_color_i, gene_mm_gs, data_subset) {
  print(module_color_i)
  ## file locations
  aracne_mod_loc = paste0("d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/",
                          data_subset, "/cytoscape_modules/aracne_", module_color_i, ".txt")
  edges_for_t20_hub_out_file = gsub(".txt$", "_top20_hub_edges.txt", aracne_mod_loc)
  if(file.exists(edges_for_t20_hub_out_file)) {
    print(paste(module_color_i, "already got edges for top 20 hubs"))
    return(module_color_i)
  }
  
  module_mm_gs = gene_mm_gs %>% 
    filter(moduleColor == module_color_i) %>% 
    select(geneSymbol, contains("GS"), contains(paste0("MM.", module_color_i, "_")))
  
  ## change names to generic
  if(data_subset == "all") {
    names(module_mm_gs) = c(names(module_mm_gs)[1:11], "Module_membership", "MM_p")
  } else if(data_subset == "ribo_w_Female") {
    names(module_mm_gs) = c(names(module_mm_gs)[1:7], "Module_membership", "MM_p")
  } else {
    names(module_mm_gs) = c(names(module_mm_gs)[1:5], "Module_membership", "MM_p")
  }  
  top_20 = module_mm_gs %>% arrange(MM_p) %>% slice(1:20) %>% select(geneSymbol) %>% unlist %>% as.character
  
  aracne_module = read_tsv(aracne_mod_loc)
  print(dim(aracne_module))
  ## from or to node in top 20 edges
  aracne_module %<>% filter((fromNode %in% top_20) | (toNode %in% top_20))
  print(dim(aracne_module))
  aracne_module %>% write_tsv(edges_for_t20_hub_out_file)
  return(module_color_i)
}

module_list = gene_mm_gs$moduleColor %>% unique
map(module_list, get_edges_for_t20_nodes, gene_mm_gs, data_subset)


##################################
## GO enrichment
##################################

expr_genes = unique(row.names(vsd))
# names(expr_genes) 
# get GO IDs from gene names
gene_map = getgo(expr_genes, 'mm9', 'ensGene')
# look for NAs in the gene_map (how many genes are represented on GO?)
sum(is.na(names(gene_map)))
sum(!is.na(names(gene_map)))
gene_map = gene_map[!is.na(names(gene_map))]


# function to return indices of genes in ARACNE module
mySelGenes = function(score) {
  return (score != 0)
}

GetTermGenesPerResult = function(go_id, gene_set_for_test, module_go) {
  # go_id = results_final[3, 1]
  term_genes = genesInTerm(module_go, go_id) %>% unlist %>% as.character
  term_genes_in_gene_set = paste(term_genes[term_genes %in% gene_set_for_test], collapse = ",")
  return(term_genes_in_gene_set)
}

GetTopTermsPerOntology = function(ontology_i, all_genes_w_scores, gene_map, gene_set_for_test, num_go_terms) {
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
  results_final = GenTable(module_go, classic = results_fet, topNodes = num_go_terms)
  results_final %<>% mutate(ontology = ontology_i)
  ## append a comma-sep list of genes in significant terms
  gene_set_term_gene_list = map(results_final$GO.ID, GetTermGenesPerResult, gene_set_for_test, module_go) %>% unlist
  results_final %<>% mutate(sig_genes = gene_set_term_gene_list)
  return(results_final)
}

### running one-offs

## take the subset of interest
gene_set_for_test = c(module_aracne_wgcna$fromNode, module_aracne_wgcna$toNode) %>% unique

# gene_set_for_test = cluster_df %>% filter(cluster_ct == 3) %>%
#   select(fromNode, toNode) %>% unlist %>% as.character %>% unique

## annotate all genes as either being or not being in the gene set of interest
all_genes_w_scores = as.numeric(expr_genes %in% gene_set_for_test)
sum(all_genes_w_scores)
length(all_genes_w_scores)
names(all_genes_w_scores) = expr_genes
# ontology_i = "MF"
num_go_terms = 5

go_results = map_df(c("CC", "MF", "BP"), GetTopTermsPerOntology, all_genes_w_scores,
                    gene_map, gene_set_for_test, num_go_terms=5)
go_results %>% select(-sig_genes)

# 
go_results %>% select(-sig_genes) %>% 
  write_tsv(paste0(wgcna_dir, "/aracne_", module_color, "_GO.txt"))


##################################
# Loop over clusters for GO
# enrichment
##################################

EnrichPerCluster = function(cluster_for_enrich, cluster_df, gene_map, expr_genes) {
  # cluster_for_enrich = 1
  print(cluster_for_enrich)
  gene_set_for_test = cluster_df %>% filter(cluster_ct == cluster_for_enrich) %>% 
    select(fromNode, toNode) %>% unlist %>% as.character %>% unique
  print(length(gene_set_for_test))
  ## annotate all genes as either being or not being in the gene set of interest
  all_genes_w_scores = as.numeric(expr_genes %in% gene_set_for_test)
  # print(sum(all_genes_w_scores))
  # print(length(all_genes_w_scores))
  names(all_genes_w_scores) = expr_genes
  go_results = map_df(c("CC", "MF", "BP"), GetTopTermsPerOntology,
                      all_genes_w_scores, gene_map, gene_set_for_test, num_go_terms=5)
  go_results %<>% mutate(cluster_ct = cluster_for_enrich)
  return(go_results)
}

cluster_go_df = map_df(1:length(clusterone), EnrichPerCluster, cluster_df, gene_map, expr_genes)

cluster_go_df %>% head %>% as.data.frame

cluster_go_df %>% 
  group_by(cluster_ct, ontology) %>% arrange(classic) %>% slice(1:2) %>% ungroup %>% 
  # select(cluster_ct, everything()) %>% head(16)
  write_tsv(paste0(wgcna_dir, "/aracne_", module_color, "_clusters_GOterms.txt"))


