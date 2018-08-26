
## D1 and D2 nuclear and whole cell RNAseq analyses

* code_pipeline
    * aracne.R: script for ARACNE co-expression calculations and hub network filtering
    * count_to_norm_fc.R and count_to_norm_functions.R: load counts data and process with TMM/voom/limma pipeline
    * de_gene_fc.R: limma DE scripts
    * deseq2_from_fc_2018_04_23.R: DESeq2 code and PCA plots
    * fastq_qc_align.py: python wrapper around alignment commands (and possibly also featurecounts if this needs heavy automating)
    * featurecounts.R: load SAM information into a counts matrix using a reference GTF or one derived from the data
    * gene_set_enrichment_gwas.R: GWAS filtering and gene set enrichment scripts
    * intron_pct_calc.R: calculate percent intron per gene
    * module_posthoc_analyses.R: plot gene kME and trait correlations
    * pca_cluster_stats.R: calculate Euclidean distance to assess relative tightness on PCA
    * variance_partition.R: calculate variance partition
    * varpart_pretty_plots.R: various custom plotting attempts for variance_partition.R output
    * wgcna.R: run WGCNA and plot results
    * wgcna_plotting_functions.R: functions for plotting WGCNA results
* gene_sets
    * human_to_mouse_orthologues.txt: ENSEMBL - gene ID maps for human and mouse genes along with confidence
    * gwas_catalog: intermediate and final files for selecting GWAS genes for enrichment tests
* README.md
    * This file
