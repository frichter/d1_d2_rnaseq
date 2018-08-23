# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 8/18/2018
# description: gene set enrichment
##############################################################

## set the home directory
setwd("D:/Dropbox/PhD/")
setwd("/Users/frichter/Dropbox (Personal)/PhD/")
setwd("/Users/felixrichter/Dropbox/PhD/")
setwd("/hpc/users/richtf01/")
options(stringsAsFactors=FALSE)

p = c("readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr")
lapply(p, require, character.only = TRUE)

## read GWAS catalog
gwa_cat = read_tsv("d1_d2_rnaseq/gene_sets/gwas_catalog/gwas_catalog_v1.0.2-associations_e93_r2018-08-14.tsv")
names(gwa_cat) = names(gwa_cat) %>% make.names

#############################################################################
# filter for interesting traits (columns: MAPPED_TRAIT, DISEASE.TRAIT)
#############################################################################

disease_grep_term = "Parkinson|Schizophrenia|Depression|Alzheimer|Bipolar|Autism"
trait_grep_term = "openness|conscientiousness|extraversion|agreeable|neurotic"
full_grep_term = paste(disease_grep_term, trait_grep_term, sep = "|")

## import additional traits to remove here (based on Hope's manual curation)
ocean_to_keep = read_tsv("d1_d2_rnaseq/gene_sets/gwas_catalog/ocean_traits_w_rm.txt") %>% 
  filter(is.na(Remove)) %>% select(-Remove)
disease_trait_to_keep = read_tsv("d1_d2_rnaseq/gene_sets/gwas_catalog/disease_traits_w_rm.txt") %>% 
  filter(is.na(Remove)) %>% select(-Remove)

trait_to_keep = bind_rows(disease_trait_to_keep, ocean_to_keep)

## only consider traits with genome-wide significance (ie P<5x10-8)
gwa_cat %<>% 
  filter(P.VALUE < 5e-8) %>% 
  inner_join(trait_to_keep)

#############################################################################
# filter for high confidence variant-gene assocotiations
#############################################################################

## Are any CNV? No
gwa_cat %>% group_by(CNV) %>% tally

# just keep non-intergenic
gwa_cat %<>% filter(INTERGENIC == 0)

## get ENSEMBL ID for mapped genes (from http://useast.ensembl.org/biomart/martview/)
## used these instructions: http://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/

ortho_tbl = read_tsv("d1_d2_rnaseq/gene_sets/human_to_mouse_orthologues.txt")
names(ortho_tbl) = names(ortho_tbl) %>% make.names
names(ortho_tbl)[5:7] = c("Mouse.orthology.confidence", "Frac_mouse_identical_to_human",
                          "Frac_human_identical_to_mouse")

## figure out the 1 to 2 human:mouse gene mappings
gwa_mouse = gwa_cat %>% 
  separate_rows(MAPPED_GENE, sep = ", ") %>% unique %>%
  mutate(lineid = 1:n()) %>%
  inner_join(ortho_tbl, by = c("MAPPED_GENE" = "Gene.name")) %>% 
  group_by(lineid) %>% mutate(human_to_mouse_mappings = n()) %>% ungroup %>% select(-lineid) %>% 
  ## for multi-mapping, only keep if Mouse.orthology.confidence == 1
  filter(!((human_to_mouse_mappings > 1) & (Mouse.orthology.confidence == 0)))

# 782/1140 mapped to mouse genes

## clean GWAS trait categories
# gwa_mouse %>% group_by(MAPPED_TRAIT, DISEASE.TRAIT) %>% tally %>% 
#   write_tsv("d1_d2_rnaseq/gene_sets/gwas_catalog/studies_per_trait.txt")

## read in simplified trait names
trait_names_simplified = read_tsv("d1_d2_rnaseq/gene_sets/gwas_catalog/trait_names_simplified.txt")
gwa_mouse %<>% inner_join(trait_names_simplified)

## how should we deal with multiple hits to the same gene? keep the best p-value per gene-trait pair
## Keep the best p-value or a meta-analysis p-value? Best p-value since w/e is reported
## could be a meta-analysis p-value. That's why Geschwind only used a few studies..
## how do multiple hits to the same gene compare (are they the same variant, for example?)
gwa_mouse %<>% 
  group_by(Trait_name, Mouse.gene.name) %>% 
  filter(P.VALUE == min(P.VALUE)) %>% 
  filter(OR.or.BETA == max(OR.or.BETA)) %>% 
  ungroup %>% 
  ## keep gwas columns you care about for downstream analyses
  ## and collapse human genes that mapped to the same mouse gene
  group_by(Mouse.gene.name, Mouse.gene.stable.ID, RISK.ALLELE.FREQUENCY,
           P.VALUE, PVALUE_MLOG, OR.or.BETA, Trait_name) %>% 
  summarise(mapped_genes_human = paste(MAPPED_GENE, collapse = ", ")) %>% ungroup %>% 
  unique

## final: 334 unique genes
gwa_mouse$Mouse.gene.name %>% unique %>% length

## export GWAS hits being used
gwa_mouse %>% select(Mouse.gene.name, Mouse.gene.stable.ID, Trait_name, P.VALUE, mapped_genes_human) %>%
  unique %>% separate_rows(Trait_name, sep = ", ") %>% unique %>% 
  group_by(Mouse.gene.name, Mouse.gene.stable.ID, Trait_name) %>% 
  # mutate(n = n()) %>% filter(n > 1)
  arrange(P.VALUE) %>% slice(1) %>% ungroup %>% dim
  # write_tsv("d1_d2_rnaseq/gene_sets/gwas_catalog/final_gwas_genes.txt")

## calculate number of traits per gene
gwa_mouse_summary = gwa_mouse %>% select(Mouse.gene.name, Mouse.gene.stable.ID, Trait_name) %>%
  unique %>% separate_rows(Trait_name, sep = ", ") %>% unique %>% 
  arrange(Mouse.gene.name, Trait_name) %>% 
  group_by(Mouse.gene.name, Mouse.gene.stable.ID) %>% 
  summarise(traits_per_gene = n(), trait_list = paste(Trait_name, collapse = ", ")) %>% 
  ungroup
  
#############################################################################
# import hub genes + neighbors, import GS_MM file
#############################################################################

## local directories:
gs_mm_loc_all = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/all/vst_bicor_signed_beta16_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"
gs_mm_loc_nuc = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/nuclear/vst_bicor_signed_beta18_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"
gs_mm_loc_wc = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/wc/vst_bicor_signed_beta14_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"
gs_mm_loc_ribo = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/ribo/vst_bicor_signed_beta9_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"
gs_mm_loc_ribo_w_F = "d1_d2_rnaseq/figures/wgcna_from_all_2018_07_12/ribo_w_Female/vst_bicor_signed_beta9_min100_mergecutheight2neg2_static99_minKMEtoStay1neg2_pamF_GS_MM_18_07_12.csv"

gs_mm_loc_list = list(gs_mm_loc_all, gs_mm_loc_nuc, gs_mm_loc_wc, gs_mm_loc_ribo, gs_mm_loc_ribo_w_F)
names(gs_mm_loc_list) = c("all", "nuclear", "wc", "ribo" , "ribo_w_Female")

## all nuclear wc ribo ribo_w_Female
data_subset = "all"
gene_mm_gs = read_csv(gs_mm_loc_list[[data_subset]])
gene_mm_gs %<>% select(-X1)

## import the hub + neighbor edge lists
hub_neighbor_file_list = list.files("d1_d2_rnaseq/manuscript/Supplementary materials/Tables/WGCNA module files/Individual files",
                                    full.names = T)
names(hub_neighbor_file_list) = hub_neighbor_file_list %>% gsub(".*/|.txt", "", .)
hub_n_nbrs = map_df(hub_neighbor_file_list, read_tsv, .id = "hub_source")

#############################################################################
# compare GWAS hits of hub genes and neighbors
#############################################################################

# overlap hub+neighbor genes with GWAS hits
# are hub genes more likely to be GWAS hits (compared to what?)
# compare hubs between methods, compare hubs to all other genes with measured RNAseq

enrich_per_trait = function(trait_i, gene_mm_gs, gwa_mouse_summary, hub_nbr_fec) {
  # print(trait_i)
  ## group by whether or not a gene is a hub+neighbor (i.e., goi) and if it's a GWAS hit
  enrich_vec = gene_mm_gs %>% 
    left_join(gwa_mouse_summary, by = c("geneSymbol" = "Mouse.gene.stable.ID")) %>% 
    mutate(Mouse.gene.name = ifelse(grepl(trait_i, trait_list), trait_list, NA)) %>% 
    # group_by(geneSymbol) %>% tally %>% filter(n > 1)
    mutate(hub_or_nbr = geneSymbol %in% hub_nbr_fec,
           gwas_hit = !is.na(Mouse.gene.name)) %>% 
    group_by(hub_or_nbr, gwas_hit) %>% tally %>% ungroup %>% 
    select(n) %>% 
    mutate(descr = c("not_goi_not_gwas", "not_goi_gwas", "goi_not_gwas", "goi_gwas")[1:n()]) %>% 
    spread(key = descr, value = n)
  ## if trait has 0 GWAS hits, add these vectors
  if(length(enrich_vec) == 2) {
    enrich_vec$not_goi_gwas = 0
  }
  if(length(enrich_vec) == 3) {
    enrich_vec$goi_gwas = 0
  }
  enrich_vec %<>% select(goi_gwas, goi_not_gwas, not_goi_gwas, not_goi_not_gwas)
  fet = fisher.test(cbind(enrich_vec[1:2] %>% t, enrich_vec[3:4] %>% t))
  enrich_vec = c(enrich_vec, "goi" = sum(enrich_vec[1:2]),
                 "fet_p" = fet$p.value, fet$estimate,
                 "ci_95_low" = fet$conf.int[[1]], 
                 "ci_95_high" = fet$conf.int[[2]]) %>% as.data.frame
  return(enrich_vec)
}

enrich_per_experiment = function(exp_source, data_subset, hub_n_nbrs, trait_list, 
                                 gene_mm_gs, gwa_mouse_summary) {
  print(exp_source)
  print(data_subset)
  ## pick the relevant GS/MM data frame
  gene_mm_gs = read_csv(gs_mm_loc_list[[data_subset]])
  gene_mm_gs %<>% select(-X1)
  
  ## subset the relevant hubs/edges of interest
  hub_nbr_fec = hub_n_nbrs %>% filter(grepl(exp_source, hub_source)) %>% select(fromNode, toNode) %>% 
    unlist %>% as.character %>% unique
  
  ## loop over al
  per_trait_enrich = map_df(trait_list, enrich_per_trait, gene_mm_gs, 
                            gwa_mouse_summary, hub_nbr_fec, .id = "trait")
  return(per_trait_enrich)
}

## get a list of GWAS traits
trait_list = gwa_mouse %>% separate_rows(Trait_name, sep = ", ") %>% select(Trait_name) %>% unique %>% 
  unlist %>% as.character
names(trait_list) = trait_list
trait_list = c("all" = "", trait_list)

exper_list = hub_n_nbrs$hub_source %>% unique
# Method_nuclear Method_RiboTag Method_whole cell
exper_list = c("Method", exper_list)
names(exper_list) = exper_list

## use the appropriate data so you have the correct background list
data_subset_list = c(rep("all", 4), rep("nuclear", 2), rep("ribo", 2), 
                     rep("ribo_w_Female", 2), rep("wc", 2))

## calculate enrichment
expr_trait_enrich = map2_df(exper_list, data_subset_list, enrich_per_experiment, 
                            hub_n_nbrs, trait_list,
                            gene_mm_gs, gwa_mouse_summary, .id = "experiment")

## write to file
# expr_trait_enrich %>% 
#   # filter(fet_p < 0.05) %>% as.data.frame
#   write_tsv("d1_d2_rnaseq/gene_sets/gwas_catalog/gwas_hit_enrichment.txt")

## plot enrichment
# https://stackoverflow.com/questions/38722202/how-do-i-change-the-number-of-decimal-places-on-axis-labels-in-ggplot2
scaleFUN = function(x) sprintf("%.2f", x)
PlotD1D2enrichment = function(experiment_i, expr_trait_enrich_d1d2) {
  print(experiment_i)
  max_ci = 25
  min_ci = 1/25
  p = expr_trait_enrich_d1d2 %>% 
    # filter(grepl(, experiment)) %>% 
    filter(experiment == experiment_i) %>%
    # filter(grepl("all", trait)) %>% 
    ## cut-off if above a certain maiximum
    mutate(ci_95_high = ifelse(ci_95_high > max_ci, max_ci, ci_95_high)) %>% 
    ## set 0 to 0.001
    mutate(ci_95_low = ifelse(ci_95_low < min_ci, min_ci, ci_95_low)) %>% 
    mutate(odds.ratio = ifelse(odds.ratio < min_ci, min_ci, odds.ratio)) %>% 
    mutate(rank_p = rank(fet_p)) %>% 
    mutate(fdr = (fet_p*rank_p)/n()) %>% 
    mutate(`FDR <0.05` = ifelse(fdr < 0.05, "Yes", "No")) %>% 
    mutate(`P <0.05` = ifelse(fet_p < 0.05, "Yes", "No")) %>% 
    ggplot(aes(x = trait, y = odds.ratio, ymax = ci_95_high, ymin = ci_95_low, col = `P <0.05`)) + 
    geom_hline(yintercept = 1, col = "grey60") +
    geom_pointrange(fatten = 1, show.legend = T, position = position_dodge(width = 0.3)) + 
    scale_color_manual(values = c("black", "red")) +
    scale_y_continuous(trans = "log2", labels = scaleFUN, limits = c(min_ci, max_ci)) +
    ylab("Odds ratio") + xlab("") +
    coord_flip() +
    theme_classic()
  p
  filename = paste0("d1_d2_rnaseq/figures/gwas_enrichment_18_08_18/gwas_enrich_", experiment_i, ".pdf")
  # print(filename)
  ggsave(filename, p, width = 4.5, height = 2.75)
}

library(qvalue)
expr_trait_enrich_d1d2 = expr_trait_enrich %>% 
  # filter(grepl("D1|D2|male", experiment)) %>% 
  # mutate(fet_p = ifelse(fet_p == 1, 0.9, fet_p)) %>% 
  # filter(fet_p == 1) %>%
  mutate(rank_p = rank(fet_p)) %>% 
  mutate(fdr = (fet_p*rank_p)/n())

# change order of x-axis for plots
level_order = expr_trait_enrich_d1d2$trait %>% unique %>% sort
level_order = rev(c(level_order[-c(5,6)], level_order[5:6]))
expr_trait_enrich_d1d2 %<>% mutate(trait = factor(trait, levels = level_order))

exp_list = expr_trait_enrich_d1d2$experiment %>% unique
# expr_trait_enrich_d1d2 %>% arrange(desc(ci_95_high))
map(exp_list, PlotD1D2enrichment, expr_trait_enrich_d1d2)

##################################
# annotate genes w GWAS results
##################################

enrich_per_experiment = function(exp_source, data_subset, hub_n_nbrs,
                                 gene_mm_gs, gwa_mouse_summary) {
  print(exp_source)
  print(data_subset)
  ## pick the relevant GS/MM data frame
  gene_mm_gs = read_csv(gs_mm_loc_list[[data_subset]])
  gene_mm_gs %<>% select(-X1)
  
  ## subset the relevant hubs/edges of interest
  hub_nbr_fec = hub_n_nbrs %>% filter(grepl(exp_source, hub_source)) %>% select(fromNode, toNode) %>% 
    unlist %>% as.character %>% unique
  gwas_mod_genes = gene_mm_gs %>% 
    select(geneSymbol, moduleColor) %>% 
    inner_join(gwa_mouse_summary, by = c("geneSymbol" = "Mouse.gene.stable.ID")) %>% 
    filter(geneSymbol %in% hub_nbr_fec)
  ## loop over al
  return(gwas_mod_genes)
}

gwas_mod_genes = map2_df(exper_list, data_subset_list, enrich_per_experiment, 
                         hub_n_nbrs, 
                         gene_mm_gs, gwa_mouse_summary, .id = "hub_neighbor_source")

gwas_mod_genes %>% 
  arrange(hub_neighbor_source, moduleColor, trait_list, Mouse.gene.name) %>% 
  write_tsv("d1_d2_rnaseq/figures/gwas_enrichment_18_08_18/hub_neighbor_gwas_intersects.txt")

gwas_mod_genes %>% 
  filter(grepl("D|Sex", hub_neighbor_source)) %>%
  arrange(Mouse.gene.name) %>% 
  group_by(Mouse.gene.name) %>% 
  # mutate(n = n()) %>% ungroup %>% filter(n >1) %>% as.data.frame
  summarise(n = n(), hub_neighbor_source = paste(hub_neighbor_source, collapse = ", ")) %>% 
  ungroup %>% 
  group_by(hub_neighbor_source) %>% tally %>% as.data.frame

### looking at the 19 genes of interest
gene_mm_gs %>% 
  left_join(gwa_mouse_summary, by = c("geneSymbol" = "Mouse.gene.stable.ID")) %>% 
  filter(geneSymbol %in% hub_nbr_fec, !is.na(Mouse.gene.name)) %>% 
  as.data.frame

# t-test of number of traits/gene between hubs vs non-hubs: NS

gene_mm_gs %>% 
  left_join(gwa_mouse_summary, by = c("geneSymbol" = "Mouse.gene.stable.ID")) %>% 
  mutate(hub_or_nbr = geneSymbol %in% hub_nbr_fec,
         gwas_hit = !is.na(Mouse.gene.name)) %>% 
  filter(gwas_hit) %>% 
  # group_by(hub_or_nbr) %>% summarise(mean(traits_per_gene), 2*sd(traits_per_gene))
  summarise(t.test(traits_per_gene ~ hub_or_nbr)$p.value)

