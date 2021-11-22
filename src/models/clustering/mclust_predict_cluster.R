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
df_sample <- df %>% sample_frac(0.01)

print("running mclust")
clust <- mclust::Mclust(df_sample %>% dplyr::select(PC1:PC25), 
                        #modelNames = "EEV",
                        G = 1:15,
                        verbose = TRUE) #5-15 used to be default, tried 1-10

saveRDS(clust, here("data/mclust/pca_mclust.Rds"))

clust_predict <- predict(clust, df %>% dplyr::select(PC1:PC25))

saveRDS(clust_predict, here("data/mclust/pca_mclust_predict.Rds"))

print("done")