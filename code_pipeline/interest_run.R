# Felix Richter, Hope Kronman
# felix.richter@icahn.mssm.edu
# 3/13/2018
# description: WGCNA
##############################################################

### Installing IntEREst on minerva
# module purge
# module load R/3.4.3
# curl -LO https://cran.r-project.org/src/contrib/DBI_1.0.0.tar.gz
# R CMD INSTALL -l . DBI_1.0.0.tar.gz 
# curl -LO https://bioconductor.org/packages/release/bioc/src/contrib/IntEREst_1.4.1.tar.gz
# R CMD INSTALL -l . IntEREst_1.4.1.tar.gz


# module purge
# module load R/3.4.3
# export TMP=/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/tmp_dir
# R

## set the home/working directory
setwd("/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq")
options(stringsAsFactors=FALSE)

library(IntEREst,lib="/hpc/users/richtf01/rLocalPackages/")
## load external libraries (order matters)
library(IntEREst)
p = c("readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr", "colorout")
lapply(p, require, character.only = TRUE)

## https://bioconductor.org/packages/release/bioc/html/IntEREst.html
## docs: https://bioconductor.org/packages/release/bioc/vignettes/IntEREst/inst/doc/IntEREst.html

gtf_file = "Mus_musculus.GRCm38.90.gtf"
testRef = referencePrepare(sourceBuild="file",
                           filePath=gtf_file, 
                           # u12IntronsChr=u12Int[,"chr"],
                           # u12IntronsBeg=u12Int[,"begin"],
                           # u12IntronsEnd=u12Int[,"end"], 
                           collapseExons=TRUE,
                           fileFormat="gtf",
                           annotateGeneIds=TRUE)

# From 6/25/2018 run:
# Warning messages:
#   1: call dbDisconnect() when finished working with a connection 
# 2: In .get_cds_IDX(type, phase) :
#   The "phase" metadata column contains non-NA values for features of type
# stop_codon. This information was ignored.

# Creating temp directory to store the results
outDir = "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/tmp_dir"
outDir = normalizePath(outDir)
# Loading suitable bam file
# possibly need to untar
# bamF = system.file(#"extdata", 
#                    "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/D1_D2_whole_cell_repeat/D1_CTRL1.sorted.bam", 
#                    package="IntEREst", mustWork=TRUE)
bamF = "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/D1_D2_whole_cell_repeat/D1_CTRL1.sorted.bam"

# possibly run getRepeatTable

# Intron retention analysis
# Reads mapping to inner introns are considered, hence 
# junctionReadsOnly is FALSE
testInterest = interest(
  bamFileYieldSize=10000,
  junctionReadsOnly=FALSE,
  bamFile=bamF,
  isPaired=TRUE,
  isPairedDuplicate=FALSE,
  isSingleReadDuplicate=NA,
  reference=testRef,
  referenceGeneNames=testRef[,"collapsed_gene_id"],
  referenceIntronExon=testRef[,"int_ex"],
  repeatsTableToFilter=c(),
  outFile=paste(outDir, "intRetRes.tsv", sep="/"),
  logFile=paste(outDir, "log.txt", sep="/"),
  # method="IntRet", ## also run with "IntSpan" or "ExEx" ## use both
  clusterNo=6,
  returnObj=TRUE,
  scaleLength= c(TRUE,FALSE, FALSE),
  scaleFragment= c(TRUE,TRUE, TRUE)
)

## 4918.371 secs for counting reads


testIntRetObj = readInterestResults(
  resultFiles= paste(outDir, "intRetRes.tsv", sep="/"), 
  sampleNames="wc_repeat_D1_CTRL1", 
  sampleAnnotation=data.frame( 
    type="wc_repeat_D1_CTRL1",
    test_ctrl="test"), 
  commonColumns=1:ncol(testRef), freqCol=ncol(testRef)+1, 
  scaledRetentionCol=ncol(testRef)+2, scaleLength=TRUE, scaleFragment=TRUE, 
  reScale=TRUE)

