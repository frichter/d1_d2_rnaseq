## Felix Richter, Hope Kronman
## 7/18/2018


## Load gene module statistics
gene_mm_gs = read_csv("d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/ribo/vst_bicor_signed_beta9_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv")
gene_mm_gs %<>% select(-X1)

## add _ to end so you can match to end of line
names(gene_mm_gs)[7:ncol(gene_mm_gs)] = paste0(names(gene_mm_gs)[7:ncol(gene_mm_gs)], "_")

## choose a module color to analyze
module_color_i = "turquoise" # "grey" # "lightgreen" # brown

## subset a specific module
module_mm_gs = gene_mm_gs %>% 
  filter(moduleColor == module_color_i) %>% 
  select(geneSymbol, contains("GS"), contains(paste0("MM.", module_color_i, "_")))

## change names to generic
names(module_mm_gs) = c(names(module_mm_gs)[1:5], "Module_membership", "MM_p")

## create -log10 p-value
module_mm_gs %<>% mutate(MM_p_neg_log10 = -log10(MM_p)) 

## plot p-value histogram
p = module_mm_gs %>% 
  ggplot(aes(x = MM_p_neg_log10)) +
  geom_histogram(bins = 40, fill = module_color_i) + 
  xlab(paste(module_color_i, "module membership\np-value (-log10)")) + ylab("Gene count") +
  theme_classic()
p
ggsave(paste0("d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/ribo/module_statistics/mm_p_hist_",
              module_color_i, ".png"), p, width = 3, height = 3)

## plot module p-value vs trait p-value 
p = ggplot(module_mm_gs, aes(x = MM_p_neg_log10, y = -log10(p.GS.D2.neurons.))) +
  geom_point(color = module_color_i) + 
  xlab("Module membership\np-value (-log10)") + 
  ylab("P-value of D2\ncorrelation(-log10)") + 
  geom_smooth(method='lm', fill = "grey80", color = "black") + 
  theme_classic()
p
ggsave(paste0("d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/ribo/module_statistics/mm_p_d2_p_", 
              module_color_i, ".png"), p, width = 3, height = 3)

## long format
gene_mm_P_long = gene_mm_gs %>% 
  select(geneSymbol, moduleColor, contains("p.MM")) %>% 
  gather(key = "MM_category", value = "MM_p", -geneSymbol, -moduleColor) %>% 
  separate(MM_category, c("t1", "t2", "MM_category"), sep = "[\\.]") %>% 
  select(-t1, -t2) %>% mutate(MM_category = gsub("_", "", MM_category)) %>% 
  filter(MM_category == moduleColor)


gene_mm_P_long %>% 
  filter(MM_p < 0.05) %>%
  group_by(moduleColor) %>% tally %>% arrange(desc(n)) %>% as.data.frame
