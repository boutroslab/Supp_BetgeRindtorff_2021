library(tidyverse)
library(here)
library(monocle3)
library(ggrastr)
library(cowplot)
library(princurve)
library(scico)
library(microbenchmark)

# I wish I could solve my path problems with the here package. 
PATH = "/dkfz/groups/shared/OE0049/B110-Isilon2/promise/"

# Loading data and formatting for easier access
obj <- readRDS(paste0(PATH, "data/interim/PhenotypeSpectrum/hdf5_umap_absolute_all_drugs.Rds"))

print("loaded data")

# defining functions
iterate_length <- function(length_in){
  print(paste0("resolution is ", length_in))
  obj_graph = cluster_cells(obj, 
                            verbose = TRUE, 
                            resolution = length_in,
                            random_seed = 1334)
  
  saveRDS(obj_graph, paste0(PATH, "data/processed/PhenotypeSpectrum/louvain_absolute_all_drugs_", length_in))
}

diagnose_l_iter <- function(length_in){
  print(paste0("diagnosing length ", length_in))
  bm <- microbenchmark(iterate_length(length_in), times = 1)
  return(bm)
}

# running benchmark
bm_l <- lapply(c(1e-8, 1e-7, 1e-6, 1e-5), diagnose_l_iter)
saveRDS(bm_l, paste0(PATH, "data/processed/PhenotypeSpectrum/benchmark_louvain_run.Rds"))