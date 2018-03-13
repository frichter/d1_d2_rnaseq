

### PCA and other QC



##############
# PCA function
##############

plotPCA = function(C, attrib, id, idx=1:2) {
  # for troubleshooting: 
  # attrib = info$library
  # id = info$Blinded.ID
  # choose which PCs to plot
  idx = sort(idx[1:2])
  # svd but only keep first idx PCs
  dcmp = svd(C, nu=max(idx), nv=0)
  # u contains the first idx eigenvectors of the expression matrix
  # d contains eigenvalues of the covariance matrix 
  percent = dcmp$d^2 / sum(dcmp$d^2)
  xlab = paste("PC", idx[1], ": ", round(percent[idx[1]]*100,1), "%", sep='')
  ylab = paste("PC", idx[2], ": ", round(percent[idx[2]]*100,1), "%", sep='')
  Covariate = droplevels(attrib)
  p_df = as.data.frame(dcmp$u[,idx])
  names(p_df) = c("x", "y")
  p_df = p_df %>% mutate(col = Covariate, ID = id)
  p = ggplot(p_df, aes(x=x, y=y, col=Covariate, label = ID)) + 
    geom_point() +
    xlab(xlab) + ylab(ylab) +
    theme_classic() #+ geom_text(data = p_df %>% filter(grepl("00425|01026|05824", ID)), 
  #            col = "black", show.legend = FALSE, check_overlap = TRUE)
  return(p)
}

vobj = readRDS("d1_d2_rnaseq/expression_data_fc/wc/vobj.RDS")
info = readRDS("d1_d2_rnaseq/expression_data_fc/wc/info.RDS")
info %<>% mutate(Cell_type = as.factor(paste(Cell_type, "neurons")))
# fit = readRDS("d1_d2_rnaseq/expression_data_fc/all/all_fit.RDS")
# R = residuals(fit, vobj)
# Mean center each gene (don't scale variance, that's what voom is for)
dataScaled = t(scale(t(vobj$E), scale = FALSE, center = TRUE))
# calculate the covariance matrix
C = cov(dataScaled)

p = plotPCA(C, info$Cell_type, info$ID, idx=1:2) ## Receptor Method
# choose to add text here:
p = p + geom_text(aes(label = ID), col = "black", show.legend = FALSE, check_overlap = F, hjust = "inward")
p = p + scale_color_manual(values =  c("green3", "red")) ## c("purple4", "skyblue")) ##
p
## neuronType method
ggsave("d1_d2_rnaseq/figures/qc_fc/voom_limma/wc_pc1_pc2_cell_type_label.png", 
       p, width = 5, height = 4, units = "in")

## repeat with residuals
R = residuals(fit, vobj)
dataScaled = t(scale(t(R), scale = FALSE, center = TRUE))
C = cov(dataScaled)


###########
# Sex check
###########

# rm_id_sex = ""
i = grep("ENSMUSG00000086503", rownames(vobj)) ## XIST
j = grep("ENSMUSG00000068457", rownames(vobj)) ## ^UTY
sum(vobj$E[j,])
vobj$E[j,] %>% head
vobj$E[i,] %>% head
sum(vobj$E[i,])#XIST = "UTY" = 
sex_gene_df = cbind(Xist = vobj$E[i,], Uty = vobj$E[j,]) %>%
  as.data.frame %>%
  mutate(file_name = colnames(vobj$E))

p = sex_gene_df %>% 
  left_join(info) %>% 
  ggplot(aes(x = Xist, y = Uty, color = gender)) + # gender , label = Sample_ID
  geom_point() +
  # geom_text(hjust = "inward", vjust = "inward") +
  # geom_text(data = sex_gene_df %>% filter(id_logic), col = "black", nudge_y = 0.3, check_overlap = TRUE) + 
  # geom_point(data = sex_gene_df %>% filter(id_logic), col = "black") +
  theme_classic()
p
ggsave("d1_d2_rnaseq/figures/qc_fc/ribo_sex_check.png", p, width = 3, height = 2.5)




