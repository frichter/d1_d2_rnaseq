# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 8/10/2017
# description: normalize count matrix RNAseq output, and integrate metadata
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

## load custom functions
# source("d1_d2_rnaseq/code_pipeline/count_to_norm_functions.R")


#############################
# FeatureCounts data cleaning
#############################

## load count matrix
fc = readRDS("d1_d2_rnaseq/expression_data_fc/counts_matrices/counts_all_17_12_06.RDS")
## gene counts from reads mapping to only exons: counts_gene_fromExons_17_12_06.RDS
## gene counts from reads mapping to introns and exons: counts_all_17_12_06.RDS

## print file names to confirm you are looking at the correct data
fc$targets

## create metadata matrix from filenames 
all_targets = fc$targets

info = as.data.frame(all_targets) %>% 
  rename(file_name = all_targets) %>% 
  ## clean the file names
  mutate(sub_name = gsub("D1_D2_", "", file_name) %>% gsub(".sam", "", .) %>% 
           gsub("CTRL", "_", .) %>% gsub("whole_cell", "wc", .) %>% 
           gsub("_repeat", "Rpt", .) %>% 
           ## add gender (male) to all previous samples
           gsub("(D[1,2])[\\.,_]", "\\1M", .) %>% 
           ## clean the ribosomal names
           gsub("([F,M])", "_\\1_", .)) %>% 
  ## create metadata columns from filename
  separate(sub_name, c("Method", "Cell_type", "gender", 
                       "sample_replicate", "Lane"), fill = "right") %>% 
  ## deal with the repeat file
  mutate(sample_replicate = ifelse(Method == "wcRpt", paste0(sample_replicate, "_2"), sample_replicate)) %>%
  mutate(Method = ifelse(Method == "wcRpt", "wc", Method)) 

## are there any counts from same sample but different lanes?
filenames_to_add = info %>% filter(!is.na(Lane)) %>% 
  group_by(Method, Cell_type, sample_replicate) %>%
  mutate(n = n()) %>% ungroup %>% filter(n > 1) %>% select(file_name) %>% unlist %>% as.character
info %<>% 
  ## deal with file with multiple lanes
  mutate(sample_replicate = ifelse(file_name %in% filenames_to_add[[2]], paste0(sample_replicate, "_L8"), sample_replicate)) %>%
  ## create a unique ID for each sample
  unite(ID, Method, Cell_type, sample_replicate, remove = F)

## how many samples are there for each?
info %>% group_by(Method, Cell_type, gender) %>% tally

#############################
# plot alignment statistics: 
# are there any outliers/abnormal samples?
#############################

stats_long = fc$stat %>% as.data.frame %>% 
  gather("file_name", "read_ct", -Status) %>% 
  ## remove lines without reads
  filter(read_ct != 0) %>% 
  ## calculate percentages
  group_by(file_name) %>% 
  mutate(read_pct = 100 * read_ct/sum(read_ct)) %>% 
  ungroup %>%  
  ## join with info
  left_join(info) 

p = stats_long %>% 
  ggplot(., aes(x = read_pct, fill = Method)) + 
  geom_histogram(bins = 25) +
  facet_wrap(~ Status, scales = "free_y") +
  scale_fill_manual(values = c("skyblue", "red3", "purple3")) + # darkblue
  theme_classic()
p
# ggsave("d1_d2_rnaseq/figures/qc_fc/mapping_features/mapping_features_pct_2018_01_29.png", 
#        p, width = 5, height = 2.5)


## check out specific files
stats_long %>% filter(grepl("Unmapped", Status)) %>% 
  # filter(file_name == "D1_D2_ribo.D1F3.sam") %>% as.data.frame
  filter(read_pct > 35)

## check out stats
stats_long %>% group_by(Status, Method) %>% 
  summarise(mean_pct = mean(read_pct), median_pct = median(read_pct), 
            mean_ct = mean(read_ct), median_ct = median(read_ct))

#############################
# subset samples here 
# if applicable
#############################

data_subset = "all" ## ribo nuclear wc all ribo_w_Female

files_to_keep = info %>% 
  # filter(grepl(gsub("_w_Female", "", data_subset), Method)) %>% ## uncomment for nuclear wc ribo
  filter(gender == "M") %>%
  ## remove the abnormal riboseq files: abnormal PCA, no DRD2, 82% of reads unassigned/unmapped
  filter(file_name != "D1_D2_ribo.D2F4.sam", 
         file_name != "D1_D2_ribo.D1F3.sam") %>% 
  select(file_name) %>% unlist %>% as.character

## check dimensions
map(fc, ~ dim(.))

fc$counts = fc$counts[, files_to_keep]
fc$targets = fc$targets[fc$targets %in% files_to_keep]
fc$stat = fc$stat[, c("Status", files_to_keep)]

## compare dimensions before and after filtering to make sure they make sense
## (ignoring counts_junction for now)
map(fc, ~ dim(.))

info %<>% filter(file_name %in% files_to_keep)

## prepare info matrix
## clean factor levels
info %<>% mutate_all(as.factor) %>% mutate_all(droplevels)
info %>% group_by(Method, Cell_type, gender) %>% tally

saveRDS(info, paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/info.RDS"))

saveRDS(fc, paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/fc_gene_from_all.RDS"))
# fc_gene_from_exon.RDS fc_gene_from_all.RDS
# fc_gene_from_exon_noF.RDS fc_gene_from_all_noF.RDS

#############################
# normalize counts with 
# voom/TMM for limma DE
#############################

## create DGEList object
x = DGEList(counts = fc$counts, genes = fc$annotation)

## filter out low-expression features (threshold may need adjustment)
## relaxed: keep if RNA cpm is greater than 1 in at least 2 samples (both abritrary)
isexpr = rowSums(cpm(x) > 1) >= 2
## stringent: keep if avg RPKM across all samples is >1
# isexpr = rowMeans(rpkm(x)) >= 1
x = x[isexpr,]

## normalize to library size within and between experiments
x = calcNormFactors(x)

## save count data
saveRDS(x, paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/norm_from_all_2018_07_23.RDS"))
## all_norm_strict.RDS all_norm_lax.RDS

## convert normalized data to RPKM and save as txt matrix
# x = readRDS(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/norm.RDS"))
# x[1:5, 1:5]
# rpkm(x)[1:5, 1:5]
write.table(rpkm(x), paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/rpkm_from_all_2018_07_23.txt"),
            quote = F, sep = "\t")

# data_subset = "all"
# rpkm_df = read.table(paste0("d1_d2_rnaseq/expression_data_fc/", data_subset, "/rpkm.txt"))

######################
# One time use
######################

## used separate counts matrices here
# cts_file_list = list.files("d1_d2_rnaseq/expression_data_fc/counts_matrices", "counts_.*.RDS", full.names = T)
# names(cts_file_list) = gsub(".*matrices/counts_", "", cts_file_list) %>% gsub("_17_11.RDS", "", .)
# cts_list = map(cts_file_list, readRDS)
# 
# 
# ## are there any counts from same sample but different lanes? If yes then could add the counts
# filenames_to_add = info %>% group_by(Method, Cell_type, sample_replicate) %>%
#   mutate(n = n()) %>% ungroup %>% filter(n > 1) %>% select(file_name) %>% unlist %>% as.character
# fc$counts[, filenames_to_add] %>% head
# new_name = gsub("L00[5-8]", "Lx", filenames_to_add) %>% unique
# new_cts = rowSums(fc$counts[, filenames_to_add])
# fc$counts %>% head
# fc$counts %<>% cbind(., new_cts)
# colnames(fc$counts)[16] = new_name


# WC highest
# filter(temp_df, (wc_rpkm > nuc_rpkm) & (wc_rpkm > ribo_rpkm))
# 
# # Higher in WC and Ribo compared to nuclear
# filter(temp_df, (wc_rpkm > nuc_rpkm) & (ribo_rpkm > nuc_rpkm))
