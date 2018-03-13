# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# description: FunPat analysis
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


# trialing FunPat
## one-time use:
# source("http://bioconductor.org/biocLite.R")
# biocLite("tseries")
# install.packages("d1_d2_rnaseq/FunPat_0.99.0.tar.gz", repos=NULL, type="source")
library(FunPat)

## look at simulated data
data("Simdata")
Simdata$replicates %>% dim
map(Simdata, head)

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

## one-to-many replicate matching
rep_1_to_many = info %>% 
  filter(Method != "wc") %>% 
  group_by(Method, Cell_type) %>% 
  mutate(main_eff = file_name[[1]]) %>% 
  slice(-1) %>% 
  ungroup %>% 
  ## arrange method so that nuclear comes before ribo
  arrange(Method) %>% 
  select(main_eff, file_name, Cell_type) %>% rename(reps = file_name)

d1_main_ids = rep_1_to_many %>% filter(Cell_type == "D1") %>% select(main_eff) %>% unique %>% unlist %>% as.character
d2_main_ids = rep_1_to_many %>% filter(Cell_type == "D2") %>% select(main_eff) %>% unique %>% unlist %>% as.character

## D1 vs D2 comparisons
case_ctrl_mapping = info %>% 
  select(Method, Cell_type, file_name) %>% 
  filter( (file_name %in% d1_main_ids) | (file_name %in% d2_main_ids) ) %>% 
  group_by(Method) %>% 
  spread(key = Cell_type, value = file_name) %>%
  ungroup %>% select(-Method)

# remove cell_type column
rep_1_to_many %<>% select(-Cell_type)

# combine all pairwise comparisons into a single DF
names(case_ctrl_mapping) = names(rep_1_to_many)
pairwise_comparison_df = bind_rows(case_ctrl_mapping, rep_1_to_many)

# combine into list (similar to Simdata)
dopa_data = list("data1" = vobj$E[, d1_main_ids], "data2" = vobj$E[, d2_main_ids], "replicates" = pairwise_comparison_df)


