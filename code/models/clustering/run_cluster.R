library(tidyverse)
library(here)
library(feather)
library(monocle3)
library(ggrastr)
library(microbenchmark)

# loading data
ods <- readRDS(here("data/monocle/bulk_pca_cluster.Rds"))
print("loaded data")

# iterating over multiple cluster sizes
#set.seed(1234)
#fraction = 0.01
#ods_sub <- ods[, sample(seq_len(ncol(ods)), ncol(ods)*fraction, replace = FALSE)]
ods_sub <- ods
#print("supsampled data")

# defining functions
iterate_k <- function(k_in){
  print(paste0("k is ", k_in))
  ods_sub_clust = cluster_cells(ods_sub, random_seed = 123, verbose = TRUE, 
                                reduction_method = "UMAP",
                                k = k_in)
  
  saveRDS(ods_sub_clust, paste0("ods_clust_", k_in))
}

diagnose_k_iter <- function(k_in){
  print(paste0("diagnosing k ", k_in))
  bm <- microbenchmark(iterate_k(k_in), times = 1)
  return(bm)
}

# running benchmark
bm_l <- lapply(c(50, 100, 1000, 10000), diagnose_k_iter)
saveRDS(bm_l, "benchmark_cluster_run.Rds")