library(tidyverse)
library(here)
library(cluster)
library(SingleCellExperiment)

# load dataset 
ods <- readRDS(here("data/monocle/bulk_pca_cluster.Rds"))

# diagnose clusters
set.seed(123456)

df <- reducedDims(ods)$PCA %>% cbind(colData(ods)) %>% as_tibble()

print("running elbow")
# elbow
# pca_kmeans_df <- tibble(k = c(1:20)) %>% 
#   mutate(kmeans = map(k, ~ kmeans(df %>% dplyr::select(PC1:PC25), centers = .x, iter.max = 10000, nstart=10, algorithm="MacQueen")),
#          var_explained = map(kmeans, ~ .x$betweenss/.x$totss))
# saveRDS(pca_kmeans_df, here("data/kmeans/pca_kmeans_df.Rds"))

# gap
print("calculating gap statistic")
gap_pca <- cluster::clusGap(df %>% dplyr::select(PC1:PC25) %>% sample_frac(0.01), FUN=kmeans, K.max=10, B=5)
saveRDS(gap_pca, here("data/kmeans/pca_gap.Rds"))

print("done")