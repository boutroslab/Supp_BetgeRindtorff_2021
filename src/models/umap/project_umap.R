# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==2) {
  stop("Arguments must be supplied (input file name, seed).n", call.=FALSE)
} 

print(args)

# seed
set.seed(args[2])

# packages
library(harmony)
library(tidyr)
library(tibble)
library(dplyr)
library(readr)
library(stringr)
library(magrittr)
library(monocle3)
library(readxl)
print("loaded libraries")

PATH = paste0(here::here(), "/")

# read object from file
obj <- readRDS(paste0(PATH, args[1]))

## Run UMAP embedding
print("starting UMAP embedding")
# harmony
obj <- reduce_dimension(obj,
                        reduction_method = "UMAP",
                        umap.min_dist = 0.1,
                        umap.n_neighbors = 15L,
                        umap.fast_sgd=TRUE,
                        cores = parallel::detectCores(),
                        verbose = TRUE)

# Save harmony result
saveRDS(obj, args[1])
print("saved UMAP projected objects at:")
print(args[1])
