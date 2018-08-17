

### deprecated code for loading data

# results_prefix = paste0("d1_d2_rnaseq/figures/wgcna_2018_04_18/", data_subset, "/vobjE_")


# alternatively use previously made voom object and fit that took cell type, method, and/or gender into account
# vobj = readRDS(paste0(home_dir, "vobj.RDS"))
# fit = readRDS(paste0(home_dir, "fit.RDS"))

# Only the residuals
# R = residuals(fit, vobj)
# datExpr = as.matrix(t(R))
# colnames(datExpr) = rownames(R)
# rownames(datExpr) = colnames(R)



# NOT residuals, but normalized (subtract each row by rowMeans, then divide each row by row SD)
# vobj_mean = rowMeans(vobj$E)
# vobj_precision_sd = 1/sqrt(vobj$weights)
# vobj_centered = vobj$E %>% as.data.frame %>% mutate_all(funs(. - vobj_mean))
# vobj_scaled = as.matrix(vobj_centered)/as.matrix(vobj_precision_sd)
# datExpr = as.matrix(t(vobj_scaled))
# colnames(datExpr) = rownames(vobj_scaled)
# rownames(datExpr) = colnames(vobj_scaled)



# wgcna_file_base = paste0( "/sc/orga/projects/chdiTrios/Felix/D1_D2_rnaseq/wgcna/onestep_",
#                           data_subset, "_signed_beta", beta_choice, "_")