# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==2) {
  stop("Arguments must be supplied (raw input file name, harmony input file name, seed).n", call.=FALSE)
} 

print(args)

# seed
set.seed(args[2])

# libraries
library(tidyr)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(here)
library(monocle3)
library(ggrastr)
library(cowplot)
library(princurve)
library(scico)
library(rhdf5) 

PATH = paste0(here::here(), "/")

# read object from file
obj <- readRDS(paste0(PATH, args[1]))

# clustering UMAP data using graph-based clustering
length_in = 1e-7
print(paste0("resolution is ", length_in))

obj = cluster_cells(obj, 
                          reduction_method = "UMAP",
                          cluster_method = "leiden",
                          verbose = TRUE, 
                          resolution = length_in,
                          random_seed = 1334)

# write object to file
obj %>% saveRDS(paste0(PATH, args[1]))