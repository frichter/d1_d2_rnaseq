


#############################
# Keep top 1k nodes
# ranked by module membership
# For modules with >1k nodes
#############################

## all nuclear wc ribo ribo_w_Female
data_subset = "ribo_w_Female"
gene_mm_gs = read_csv(gs_mm_loc_list[[data_subset]])

gene_mm_gs %<>% select(-X1)

large_module_list = gene_mm_gs %>% group_by(moduleColor) %>% tally %>% 
  filter(n > 1000) %>% select(moduleColor) %>% unlist %>% as.character

large_module_list

## add _ to end so you can match to end of line
names(gene_mm_gs)[7:ncol(gene_mm_gs)] = paste0(names(gene_mm_gs)[7:ncol(gene_mm_gs)], "_")

module_color_i = large_module_list[[1]] # "turquoise" # "grey" # "lightgreen" # brown

filter_top1k_nodes = function(module_color_i, gene_mm_gs, data_subset) {
  print(module_color_i)
  
  module_mm_gs = gene_mm_gs %>% 
    filter(moduleColor == module_color_i) %>% 
    select(geneSymbol, contains("GS"), contains(paste0("MM.", module_color_i, "_")))
  
  ## change names to generic
  if(data_subset == "all") {
    names(module_mm_gs) = c(names(module_mm_gs)[1:11], "Module_membership", "MM_p")
  } else {
    names(module_mm_gs) = c(names(module_mm_gs)[1:5], "Module_membership", "MM_p")
  }  
  top_1k = module_mm_gs %>% arrange(MM_p) %>% slice(1:1000) %>% select(geneSymbol) %>% unlist %>% as.character
  
  aracne_mod_loc = paste0("d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/",
                          data_subset, "/cytoscape_modules/aracne_", module_color_i, ".txt")
  
  aracne_module = read_tsv(aracne_mod_loc)
  print(dim(aracne_module))
  aracne_module %<>% filter(fromNode %in% top_1k, toNode %in% top_1k)
  print(dim(aracne_module))
  aracne_module %>% write_tsv(gsub(".txt$", "_top1k.txt", aracne_mod_loc))
  return(module_color_i)
}

map(large_module_list, filter_top1k_nodes, gene_mm_gs, data_subset)

###########################
# Figuring out radiality
###########################

# https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-8-S4-S11


###########################
# Figuring out radiality
###########################


# module_aracne_wgcna %>% 
#   select(-mic_aracne, -direction) %>% 
#   write_tsv(paste0(wgcna_dir, "/aracne_", module_color, "_wgcna_weights.txt"), col_names = F)


# module_aracne_wgcna %>% 
#   select(-weight_wgcna, -mic_aracne, -direction) %>% 
#   write_tsv(paste0(wgcna_dir, "/aracne_", module_color, "_noWeights.txt"))

cONE_cmd = paste0("time java -jar /mnt/c/Program\\ Files/cluster_one-1.0.jar -F plain ",
                  "/mnt/d/Dropbox/PhD/", wgcna_dir, "/aracne_", module_color, "_wgcna_weights.txt > ",
                  "/mnt/d/Dropbox/PhD/", wgcna_dir, "/aracne_", module_color, "_clusterOne_rawOutput.txt")
write(cONE_cmd, paste0(wgcna_dir, "/../clusterOne_cmd.txt"), append=T)
## running clusterONE
# cd /mnt/d/Dropbox/PhD/d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/nuclear/cytoscape_modules
## java -jar /mnt/c/Program\ Files/cluster_one-1.0.jar --version
# java -jar /mnt/c/Program\ Files/cluster_one-1.0.jar 
# time java -jar /mnt/c/Program\ Files/cluster_one-1.0.jar -F \
# plain aracne_turquoise_wgcna_Weights.txt > aracne_turquoise_wgcna_Weights_clusterone.txt
# cut -f1,2,5 aracne_fromR_turquoise_t50.txt | sed 1d > aracne_turquoise_wgcna_t50.txt
# aracne_turquoise_wgcna_t50_clusterone.txt

module_aracne_wgcna = read_tsv(paste0(wgcna_dir, "/aracne_", module_color, ".txt"))

############## FIGURE OUT THE MINIMUM DENSITY FLAGS, need --min-density 0.25 for ribo brown

###########################
# after running clusterOne,
# keep the subset of edges
# that form tight networks
###########################

# read as list
clusterone = readLines(paste0(wgcna_dir, "/aracne_", module_color, "_clusterone_rawOutput.txt"))

# cluster_ct = 1

getClusterWeights = function(cluster_ct, clusterone) {
  cluster_i = clusterone[[cluster_ct]]
  # find the edges that correspond to a cluster
  cluster_vec_i = strsplit(cluster_i, "\t")[[1]]
  
  cluster_i_weights = module_aracne_wgcna %>% 
    filter(fromNode %in% cluster_vec_i,
           toNode %in% cluster_vec_i) %>% 
    mutate(cluster_ct = cluster_ct)
  print("Observed/max number of edges (confirm <=1)")
  print(nrow(cluster_i_weights)/ (length(cluster_vec_i)^2))
  return(cluster_i_weights)
}

cluster_df = map_df(1:length(clusterone), getClusterWeights, clusterone)


drd1 = "ENSMUSG00000021478"
drd2 = "ENSMUSG00000032259"

## Genes can belong to multiple clusters
cluster_df %>% 
  # filter((fromNode %in% drd1) | (toNode %in% drd1)) %>%
  # filter((fromNode %in% drd2) | (toNode %in% drd2)) %>% 
  group_by(cluster_ct) %>%
  summarise(uniq_from = unique(fromNode) %>% length, uniq_to = unique(toNode) %>% length, n()) %>% 
  ungroup %>% 
  summarise(sum(uniq_from), sum(uniq_to))


cluster_df %>% write_tsv(paste0(wgcna_dir, "/aracne_", module_color, "_clusters.txt"))


###########################
# load aracne/wgcna module
# to decrease size
###########################

# module_aracne_wgcna = read_tsv(paste0(wgcna_dir, "/aracne_fromR_", module_color, ".txt"))
# 
# # any other filters to decrease module size?
# module_aracne_wgcna %>% filter(!is.na(weight_wgcna)) %>% dim
# # ggplot(module_aracne_wgcna, aes(x = mic_aracne)) + geom_density() + theme_classic()
# # ggplot(module_aracne_wgcna, aes(x = weight_wgcna)) + geom_density() + theme_classic()
# # ggplot(module_aracne_wgcna, aes(x = weight_wgcna, y = mic_aracne)) + geom_point() + theme_classic()
# 
# # remove the bottom quartile of edges
# module_aracne_wgcna_top = module_aracne_wgcna %>% 
#   filter(weight_wgcna > quantile(weight_wgcna, 0.5))
#          # mic_aracne > quantile(mic_aracne, 0.25))
# 
# module_aracne_wgcna %>% #head(100000) %>% 
#   select(weight_wgcna, mic_aracne) %>% cor
# 
# # write to file (without weights to minimize size)
# module_aracne_wgcna_top %>% 
#   select(-weight_wgcna, -mic_aracne, -direction) %>% 
#   write_tsv(paste0(wgcna_dir, "/aracne_fromR_", module_color, "_wgcna.txt"))
