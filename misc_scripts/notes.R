

### Notes


## 8/10/2017
## downloaded the gencode GTFs here
## http://www.gencodegenes.org/mouse_releases/

source("https://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
# http://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html
# check related packages for enrichment

library(org.Mm.eg.db)
library(annotate)

geneIDs = as.character(row.names(vobj$E[i, ]))
geneSyms = unlist(lookUp(t30$gene_name, 'org.Mm.eg', 'ENSEMBL'))
mapIds(org.Mm.eg.db, t30$gene_name, "SYMBOL","ENSEMBL")

## 2/21/2018

# trialing FunPat
source("http://bioconductor.org/biocLite.R")
biocLite("tseries")
install.packages("d1_d2_rnaseq/FunPat_0.99.0.tar.gz", repos=NULL, type="source")
library(FunPat)



SEL.TS.AREA

### set file names
data_subset = "all" ## all nuclear wc ribo
home_dir = paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/")
x_loc = paste0(home_dir, "norm.RDS")
info_loc = paste0(home_dir, "info.RDS")
vobj_loc = paste0(home_dir, "vobj.RDS")

## read normalized count and metadata matrices
x = readRDS(x_loc) 
info = readRDS(info_loc) 
vobj = readRDS(vobj_loc) 

vobj$E %>% head
x$counts %>% head

info %>% 
  filter(Method != "wc") %>% 
  group_by(Method, Cell_type) %>% 
  mutate(main_eff = file_name[[1]]) %>% 
  slice(-1) %>% 
  ungroup %>% 
  select(main_eff, file_name)

nuc_d1 = info %>% filter(Method == "nuclear", Cell_type == "D1") %>% 
  slice(1) %>% select(file_name) %>% unlist %>% as.character
nuc_d1_replicates = info %>% filter(Method == "nuclear", Cell_type == "D1") %>% 
  slice(-1) %>% select(file_name) %>% unlist %>% as.character

info %>% 
  group_by

ribo_A = info %>% filter(method == "ribo")
ribo_replicates = info %>% filter(method == "ribo")

replicates_tbl = as.data.frame()




