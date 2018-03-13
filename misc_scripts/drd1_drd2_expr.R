setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

## load external libraries (order matters)
p = c("limma", 
      ## generic dataprocessing packages:
      "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

## looking at D1 and D2 expression across all datasets
drd1 = "ENSMUSG00000021478"
drd2 = "ENSMUSG00000032259"

# vobj = readRDS("d1_d2_rnaseq/expression_data/all_vobj.RDS")
# info = readRDS("d1_d2_rnaseq/expression_data/all_info.RDS")
d_expr = vobj$E[c(drd1, drd2), ] %>% as.data.frame %>% 
  mutate(ens_id = c(drd1, drd2)) %>% 
  gather(key = "file_name", value = "expr", -ens_id) %>% ## ID
  mutate(gene = ifelse(ens_id == drd1, "DRD1", "DRD2")) %>% 
  left_join(info) %>% 
  # mutate(Receptor = paste(Receptor, "neurons"))
  mutate(Cell_Type = paste(Cell_type, "neurons"))

p = ggplot(d_expr, aes(x = Method, y = expr, col = gene)) +
  geom_point() + 
  facet_wrap(~ Cell_Type) +
  scale_color_manual(values = c("green3", "red")) + 
  ylab("Expression (log2 cpm)") +
  theme_classic()
p
ggsave("d1_d2_rnaseq/figures/individual_genes/drd1_drd2.png", p, width = 4, height = 3)
