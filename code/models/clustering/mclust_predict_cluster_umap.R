library(tidyverse)
library(here)
library(cluster)
library(SingleCellExperiment)
library(mclust)

# load dataset 
ods <- readRDS(here("vignettes/07_organoid_unsupervised_exploration/monocle/bulk_pca_cluster.Rds")) #used to be bulk_pca_cluster_strat.Rds

# diagnose clusters
set.seed(123456)

df <- reducedDims(ods)$UMAP %>% cbind(colData(ods)) %>% as_tibble()
df_sample <- df %>% group_by(Line, Replicate) %>% sample_frac(0.1) %>% ungroup()

print("running mclust")
clust <- mclust::Mclust(df_sample %>% dplyr::select(V1:V2), 
                        #modelNames = "EVV", 
                        G = 1:25,
                        verbose = TRUE) #5-15 used to be default, tried 1-10

saveRDS(clust, here("vignettes/07_organoid_unsupervised_exploration/mclust/umap_mclust.Rds"))

clust_predict <- predict(clust, df %>% dplyr::select(V1:V2))

#saveRDS(clust_predict, here("data/mclust/umap_mclust_predict)
#saveRDS(clust_predict, here("data/mclust/umap_mclust_predict_k6.Rds"))
saveRDS(clust_predict, here("vignettes/07_organoid_unsupervised_exploration/mclust/umap_mclust_predict.Rds"))

print("done")