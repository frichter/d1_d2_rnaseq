

PrepNucMetadata = function(nuclear_cts) {
  info_nuc = as.data.frame(names(nuclear_cts)[-1])
  names(info_nuc) = "ID"
  info_nuc %<>% 
    ## extract covariates from names
    separate(ID, into = c("Receptor", "Replicate"), sep = "\\.", remove = F) %>% 
    ## convert to factors (change to mutate_at if adding continuous variables)
    mutate_all(as.factor) %>% mutate_all(droplevels)
  return(info_nuc)
}


PrepWCMetadata = function(whole_cell_cts) {
  info_wc = as.data.frame(names(whole_cell_cts)[-1])
  names(info_wc) = "ID"
  info_wc %<>% 
    ## extract covariates from names
    separate(ID, into = c("Receptor", "Replicate"), sep = "CTRL", remove = F) %>% 
    ## convert to factors
    mutate_all(as.factor) %>% mutate_all(droplevels)
  return(info_wc)
}
