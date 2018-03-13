# Felix Richter
# felix.richter@icahn.mssm.edu
# 10/31/2017
# Description: count number of features (including introns) in SAM files
################################################################################

## Minerva modules
# module load R/3.4.1
# R

## load R libraries:
p = c("Rsubread", "colorout", "magrittr", "purrr", "stringi", "dplyr",
      "ggplot2", "tidyr", "readr")
lapply(p, require, character.only = TRUE)

options(stringsAsFactors=FALSE)

## go to working directory
setwd("/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq")

## testing
# input_file = "D1_D2_odd_repeat/D1_CTRL1.sam"
# stringtie_gtf = "D1_D2_odd_repeat/D1_CTRL1_refguided.gtf"

## directories
## D1_D2_whole_cell D1_D2_whole_cell_repeat D1_D2_nuclear
input_dir_list = c("D1_D2_whole_cell", "D1_D2_whole_cell_repeat", "D1_D2_nuclear", "D1_D2_ribo")
input_file_list = map(input_dir_list, ~ list.files(., ".*sam", full.names = T)) %>% unlist

## actual featureCounts command
fc = featureCounts(files = input_file_list,
                   # annot.inbuilt = "mm10",
                   annot.ext = "Mus_musculus.GRCm38.90.gtf",
                   ## options: "Mus_musculus.GRCm38.90.gtf" stringtie_gtf
                   isGTFAnnotationFile = T,
                   ## if using introns, use "gene" for GTF.featureType
                   ## for no introns, use "exon"
                   GTF.featureType = "exon", ## Checked GTF that this has introns
                   GTF.attrType = "gene_id",
                   juncCounts = T,  ## provides a fc$counts_junction matrix
                   useMetaFeatures=FALSE, ## TRUE for gene-level, FALSE for exons
                   allowMultiOverlap = TRUE,
                   nthreads = 12,
                   strandSpecific = 0, #0, 1 or 2
                   countMultiMappingReads = FALSE,
                   isPairedEnd = TRUE,
                   tmpDir = "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/tmp_dir",
                   verbose = TRUE)


## 20% when using GTF.featureType = "exon", 46% with "gene"
## when using Mus_musculus.GRCm38.90.gtf
## over half the reads are intronic

## counts_WC_17_11 counts_Nuc_17_11
saveRDS(fc, file = "counts_matrices/counts_exon_17_12_08.RDS")


## manuals
## http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf
## http://bioconductor.org/packages/release/bioc/manuals/Rsubread/man/Rsubread.pdf
