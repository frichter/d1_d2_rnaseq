
## D1 and D2 nuclear and whole cell RNAseq analyses

* code_pipeline
    * _count_to_norm.R and count_to_norm_functions.R:_
        * load counts data: either HTSeq, which used an unknown GTF, or FeatureCounts, which used gencode.vM13.annotation.gtf.gz
    * _de_gene.R_
    * _featurecounts.R:_ load SAM information into a counts matrix using a reference GTF or one derived from the data
    * _pca_sex_check.R:_ standard QC
    * _qc_align_count.py:_ python wrapper around alignment commands (and possibly also featurecounts if this needs heavy automating)
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
    * _drd1_drd2_expr.R:_ checks normalized expression of D1 and D2 receptor genes
    * _notes.R, notes.sh:_ log of activity in R and shell scripts
* README.md
    * This file.. :)
* reference_data
    * _gencode.vM13.annotation.gtf.gz:_ mouse reference annotation with multiple feature layers (e.g., gene, exon)
