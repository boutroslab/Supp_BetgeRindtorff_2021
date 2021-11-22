library(tidyverse)
library(here)
#library(feather)
library(monocle3)
#library(ggrastr)
library(microbenchmark)

# loading data
ods_sub <- readRDS(here("data/monocle/bulk_pca_cluster_subsampled.Rds"))
print("loaded clustered data")

# accidently introduced
#k_in = 100
ods_sub = cluster_cells(ods_sub, random_seed = 123, verbose = TRUE, 
                        reduction_method = "UMAP")

# defining functions
iterate_length <- function(length_in){
  print(paste0("k is ", length_in))
  ods_sub_graph = learn_graph(ods_sub, verbose = TRUE, 
                        close_loop = TRUE, 
                        learn_graph_control = list(minimal_branch_len = length_in))
  
  saveRDS(ods_sub_graph, paste0("ods_sub_graph_", length_in))
}

diagnose_l_iter <- function(length_in){
  print(paste0("diagnosing length ", length_in))
  bm <- microbenchmark(iterate_length(length_in), times = 1)
  return(bm)
}

# running benchmark
bm_l <- lapply(c(10, 20, 30, 40, 50, 100, 1000, 10000), diagnose_l_iter)
saveRDS(bm_l, "benchmark_graph.Rds")