
# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 6/10/2018
# description: Plot distribution of module weights
##############################################################

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("limma", "edgeR", "annotate", "org.Mm.eg.db", "DESeq2",
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)


module_list = list.files("d1_d2_rnaseq/figures/wgcna_2018_04_18/ribo/cytoscape_modules",
                         "cs_edges.*", full.names = T)
grep("blue", module_list)
wgcna_network = read_tsv(module_list[[2]])


# p = ggplot(wgcna_network, aes(x = weight)) +
#   geom_density() +
#   theme_classic()
# p

##

net_wide = wgcna_network %>% 
  select(-fromAltName, -toAltName, -direction) %>% 
  spread(key = toNode, value = weight, fill = 0)

net_wide[1:5, 1:5]
sum(is.na(net_wide))

net_mtx = net_wide %>% select(-fromNode) %>% as.matrix
row.names(net_mtx) = net_wide$fromNode
id_order = match(rownames(net_mtx), colnames(net_mtx))
rownames(net_mtx)[!(rownames(net_mtx) %in% colnames(net_mtx))]
colnames(net_mtx)[!(colnames(net_mtx) %in% rownames(net_mtx))]
net_mtx[]
net_mtx[diag(net_mtx) != 0, diag(net_mtx) != 0]

net_wide
aracne_out_test = aracne(net_wide)



