library(tidyverse)
library(here)
library(cluster)
library(SingleCellExperiment)
library(mclust)

# load dataset 
ods <- readRDS(here("data/monocle/bulk_pca_cluster.Rds"))

# diagnose clusters
set.seed(123456)

df <- reducedDims(ods)$PCA %>% cbind(colData(ods)) %>% as_tibble()

print("running mclust")
clust <- mclust::Mclust(df %>% dplyr::select(PC1:PC25), 
                        #modelNames = "EEV",
                        G = 1:15,
                        verbose = TRUE) #5-15 used to be default, tried 1-10

saveRDS(clust, here("data/mclust/pca_mclust.Rds"))

print("running mclust on umap")

df <- reducedDims(ods)$UMAP %>% cbind(colData(ods)) %>% as_tibble()

clust <- mclust::Mclust(df %>% dplyr::select(V1:V2), 
                        #modelNames = "EEV",
                        G = 1:15,
                        verbose = TRUE) #5-15 used to be default, tried 1-10

saveRDS(clust, here("data/mclust/umap_mclust.Rds"))
print("done")