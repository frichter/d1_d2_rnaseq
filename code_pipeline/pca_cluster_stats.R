

## continuation of deseq2_from_fc_2018_04_23.R (start there first to generate p_data)

## centroid clustering

p_data %>% head

## calculate centroids per cluster
p_centroids = p_data %>% group_by(group) %>% 
  mutate(centroid_pc1 = mean(PC1), centroid_pc2 = mean(PC2)) %>% 
  ungroup

## calculate average Euclidean distance to centroid
p_centroids %<>% 
  rowwise() %>% 
  mutate(dist_euclid = sqrt((PC1 - centroid_pc1)^2 + (PC2 - centroid_pc2)^2)) %>% 
  ungroup

## check work (should equal 3.73)
dist(rbind(c(-28.3, -12.4), c(-31.3, -10.2)))

## graph and calculate significance of difference between centroids
p = p_centroids %>% 
  group_by(Method) %>% 
  mutate(mean_dist = mean(dist_euclid)) %>% ungroup %>% 
  ggplot(aes(x = dist_euclid), color = "black") +
  geom_histogram(bins = 20) +
  facet_wrap(~Method, ncol = 1) +
  scale_y_continuous(breaks = c(0, 5, 10)) +
  geom_vline(aes(xintercept = mean_dist), color = "red") +
  xlab("PC1/PC2 Euclidean distance\nto cluster centroid") +
  theme_classic()
p
filename = paste0("d1_d2_rnaseq/figures/qc_fc/deseq2/", data_subset, 
                  "_from_", data_source, "_centroid_distance_2018_08_09.svg")
ggsave(filename, p, width = 2.5, height = 3, units = "in")

## check statistical significance of difference
fit = aov(dist_euclid ~ Method, p_centroids %>% mutate(Method = factor(Method)))
summary(fit)
TukeyHSD(fit)

## ribo vs nuclear
t.test(dist_euclid ~ Method, p_centroids %>% filter(Method !="wc"))
## ribo vs wc
t.test(dist_euclid ~ Method, p_centroids %>% filter(Method !="nuclear"))
## nuclear vs wc
t.test(dist_euclid ~ Method, p_centroids %>% filter(Method !="ribo"))


p_centroids %>% 
  group_by(Method) %>% 
  summarise(mean_dist = mean(dist_euclid))

## calculate distance between cluster centroids
p_method_centroids = p_centroids %>% select(contains("centroid")) %>% unique %>% as.matrix
row.names(p_method_centroids) = unique(p_centroids$Method)
dist(p_method_centroids)

## Mahalanobis distance: distance between a point and distribution

## plot centroids on PCA
p = p_centroids %>%
  ggplot(., aes(PC1, PC2, color=group)) +
  geom_point(size=2.5) + ## 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  # coord_fixed() +
  ### order is: nuclear, ribo, wc
  scale_color_manual(values = c("royalblue", "red", "black")) + ## c("skyblue", "red3", "purple3")
  # scale_color_manual(values =  c("green3", "red")) + ##c("blue", "grey")) + ##
  # geom_text(aes(label = group), col = "black", show.legend = FALSE, check_overlap = F, hjust = "inward") +
  geom_point(mapping = aes(x=centroid_pc1, y=centroid_pc2), color = "yellow3", size = 3, alpha=0.5) +
  theme_classic()
p
filename = paste0("d1_d2_rnaseq/figures/qc_fc/deseq2/", data_subset, 
                  "_from_", data_source, "_pc1_pc2_method_w_centroid_2018_08_09.svg")
## _pc1_pc2_gender.png _pc1_pc2_cell_type.png _pc1_pc2_method_2018_06_23.png
ggsave(filename, p, width = 3.5, height = 2.5, units = "in")
