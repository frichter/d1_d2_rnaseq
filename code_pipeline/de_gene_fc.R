# Felix Richter, Hope Kronman
# Created 8/10/2017
# Description: run DE analysis with voom/limma
##############################################################

setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("limma", "edgeR", 
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

### set file names
data_subset = "nuclear" ## ribo nuclear wc all ribo_w_Female
data_source = "all" ## exon all
home_dir = paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/")
x_loc = paste0(home_dir, "norm_from_", data_source, ".RDS") ## _norm_lax.RDS _norm_strict.RDS
info_loc = paste0(home_dir, "info.RDS")
vobj_loc = paste0(home_dir, "vobj_", data_source, ".RDS")
fit_loc = paste0(home_dir, "fit_", data_source, ".RDS")

## read normalized count and metadata matrices
x = readRDS(x_loc) 
info = readRDS(info_loc) 

## looking at the chromosomes that were captured:
x$genes %>% select(Chr) %>% as.data.frame %>% unlist %>% as.character %>%
  paste(collapse = ";") %>% strsplit(., ";") %>% unlist %>% unique
# x$genes %>% filter(grepl("MT", Chr)) %>% tally
## Looks like non-nuclear DNA https://www.ncbi.nlm.nih.gov/nuccore/JH584304.1

## create model matrix
design = model.matrix(~ Cell_type + 0, info) ## Cell_type + Method + gender

## estimate variance as function of mean expression
vobj = voom(x, design, plot=TRUE)
# saveRDS(vobj, vobj_loc)

## update vobj to subtract out the median and IQR per sample

col_meds = colMedians(vobj$E)
col_iqrs = colIQRs(vobj$E)

expr_norm = apply(vobj$E, 1, function(x) (x - col_meds)/col_iqrs)
vobj$E = t(expr_norm)
rownames(vobj$E) = colnames(expr_norm)
colnames(vobj$E) = rownames(expr_norm)

## fit a linear model to each gene, using the design matrix
fit = lmFit(vobj, design)
# saveRDS(fit, fit_loc)

#### simple DE (comment out lines you don't want to compare)
cont_matrix = makeContrasts(Cell_type_effect = Cell_typeD2 - Cell_typeD1, 
                            # Method_effect = Methodwc, ##
                            levels = design)
# cont_matrix = makeContrasts(levels = design)
fit2 = contrasts.fit(fit, cont_matrix)
fit2 = eBayes(fit2)

summary(decideTests(fit2, p.value = 0.05, lfc = 0.38))

## Cell_type_effect Receptor_effect Method_effect
topSet = topTable(fit2, coef = "Cell_type_effect", p.value = 1, number = nrow(fit2))

### sanity checks that D1 and D2 are differentially expressed
drd1 = "ENSMUSG00000021478"
drd2 = "ENSMUSG00000032259"

topSet[drd1, ]
topSet[drd2, ] ## nuclear D2 is not differentially expressed fit2[drd2, ] 

# topSet %>% as.data.frame %>% 
#   mutate(gene_id = row.names(topSet)) %>% 
#   write_tsv(., paste0("d1_d2_rnaseq/de_tables/fc_voom_limma/", data_subset, "_d1_v_d2_2018_06_17.txt"))





######################
# One time use
######################

# sum(!(row.names(x_nuc) %in% row.names(x_wc)))
# sum(!(row.names(x_wc) %in% row.names(x_nuc)))

all_cor = cor(vobj$E)

min(all_cor)



