library(tidyverse)
library(here)
library(monocle3)

set.seed(21231)
PATH = "/dkfz/groups/shared/OE0049/B110-Isilon2/promise/"

h_line = c("D004T01", "D007T01", "D010T01", "D013T01", "D018T01", "D020T01", "D020T02", "D021T01", "D022T01", "D027T01", "D030T01", "D046T01", "D054T01", "D055T01", "D019T01")

for(loi in h_line){
  plane <- read_csv(paste0(PATH, "data/interim/FeatureAnalysis/drug_effects/human/",loi,"/SVM_Profiles_mean_PCA_", loi, "_25components.csv")) 
  print("read lm vectors")
  plane <- plane %>% 
    mutate(id = paste0(loi, "_", X1)) %>%
    dplyr::select(-X1) %>% 
    as.data.frame() %>% 
    column_to_rownames("id") %>% 
    as.matrix()
  
  activity <- read_csv(paste0(PATH, "data/interim/FeatureAnalysis/drug_effects/human/",loi,"/SVM_Accuracies_PCA_",loi, "_25components.csv")) %>% 
    janitor::clean_names() %>% 
    dplyr::rename(drug = x1) %>% 
    mutate(line = loi,
           id = paste0(line, "_", drug)) %>% 
    as.data.frame() %>% 
    column_to_rownames("id")
  
  index = which(activity$auc_mean >= 0.85)
  
  ## combining annotation data with harmonized PCA information
  obj <- new_cell_data_set(plane[index,] %>% t(),
                           cell_metadata = activity[index,])
  ## I manually inject the PCA compression of the data into the object
  reducedDims(obj)$PCA <- plane[index,]
  
  print("running gridsearch")
  for(i_dist in c(0.05, 0.1, 0.5, 1))
  {
    for(j_nn in c(5, 15, 30))
    {
      for(k_res in c(1e-2, 1e-1, 1)){
        print(i_dist)
        print(j_nn)
        print(k_res)
        ## Run UMAP embedding
        obj <- reduce_dimension(obj,
                                reduction_method = "UMAP",
                                umap.min_dist = i_dist,
                                umap.n_neighbors = j_nn,
                                umap.fast_sgd=TRUE, 
                                cores = parallel::detectCores(),
                                verbose = TRUE)
        obj = cluster_cells(obj, resolution=k_res)
        # Save overview plot
        plot_cells(obj) + ggsave(paste0(PATH, "reports/figures/lm_umap", loi, "_", i_dist, "_", j_nn, "_", k_res,".pdf"))
        # Save harmony result
        saveRDS(obj, paste0(PATH, "data/processed/PhenotypeSpectrum/drug_effects/lm_umap", loi, "_", i_dist, "_", j_nn, "_", k_res, ".Rds"))
      }
    }
  }
}
