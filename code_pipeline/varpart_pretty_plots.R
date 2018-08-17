# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 1/29/2018
# description: Variance partition pretty plotting
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


varPart = readRDS("d1_d2_rnaseq/expression_data_fc/all/varpart_from_all_2018_07_23.RDS")
varPart = sortCols(varPart)

colnames(varPart) = c("Method", "Cell type", "Residuals") # "Gender",

varPart %<>% 
  mutate(gene_names = row.names(varPart)) %>% 
  gather(key = "covar", value = "var_explained", -gene_names) %>% 
  mutate(var_explained = var_explained*100) %>% 
  mutate(covar = factor(covar, levels = colnames(varPart)))

RColorBrewer::display.brewer.pal(8, "Dark2")

varpart_colors = RColorBrewer::brewer.pal(8, "Dark2")[c(3,6,8)]
varpart_colors
varpart_colors_long = plyr::mapvalues(varPart$covar, unique(varPart$covar), varpart_colors)
## 10% darker found here under "SHADES OF #": https://www.colorbook.io/hexcolors/view/666666
varpart_colors_darker = c("#6C539D", "#B38502", "#4C4C4C")
p = varPart %>% 
  ggplot(aes(x = covar, y = var_explained)) +
  geom_violin(aes(fill = covar, color = covar), scale = "width", trim = T, adjust = 0.75) +
  geom_boxplot(aes(color = covar), width = 0.2, show.legend = F, outlier.size = 0.35) +
  scale_fill_manual(values = varpart_colors) + 
  scale_color_manual(values = varpart_colors_darker) + 
  xlab("") + ylab("Variance Explained (%)") +
  theme_classic()
p
ggsave("d1_d2_rnaseq/figures/variance_partition_fc_2018_07_23/var_explained_darks2_box_small_outs.svg",
       p, width = 3.75, height = 3.25)
#  variance_partition_fc_2018_06_23

######## Testing breaking up into quantils
# color_scale = colorpanel(100, "Blue", "Black", "Yellow")[c(1:25, 45:55, 75:100)]
# color_scale = c(rep(color_scale[[1]], 50), color_scale, rep(color_scale[[length(color_scale)]], 60))
# 
# p = varPart %>% 
#   # group_by(covar) %>% summarise(mean = mean(var_explained))
#   ggplot(aes(x = covar, y = var_explained, fill = covar)) + # , color = var_explained
#   geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
#   # geom_boxplot(width = 0.2, show.legend = F) +
#   # geom_point() +
#   scale_fill_manual(values = varpart_colors) +
#   xlab("") + ylab("Variance Explained (%)") +
#   theme_classic()
# p
# 
# ggsave("d1_d2_rnaseq/figures/variance_partition_fc_2018_06_23/var_explained_col_by_y.png",
#        p, width = 3.75, height = 3.25)

