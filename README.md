
## D1 and D2 nuclear and whole cell RNAseq analyses

* code_pipeline
    * fastq_qc_align.py: python wrapper around alignment commands (and possibly also featurecounts if this needs heavy automating)
    * count_to_norm.R and count_to_norm_functions.R
        * load counts data: either HTSeq, which used an unknown GTF, or FeatureCounts, which used gencode.vM13.annotation.gtf.gz
    * de_gene.R
    * featurecounts.R: load SAM information into a counts matrix using a reference GTF or one derived from the data
    * pca_sex_check.R: standard QC
* de_tables
* expression_data_fc
    * input, intermediate, and processed data for _FeatureCounts_
    * counts_matrices: count matrices
    * normalized_matrices_stringent: normalized matrix after filtering for mean(RPKM(gene)) >= 1
* expression_data_htseq
    * input, intermediate, and processed data for _HTSeq counts_
* figures
* misc_scripts
    * One time or small analysis scripts that are not part of main pipelines. Might be incorporated into pipelines if they see heavy use
    * drd1_drd2_expr.R: checks normalized expression of D1 and D2 receptor genes
    * notes.R, notes.sh: log of activity in R and shell scripts
* README.md
    * This file.. :)
* reference_data
    * _gencode.vM13.annotation.gtf.gz:_ mouse reference annotation with multiple feature layers (e.g., gene, exon)
