#!/usr/bin/Rscript

# packages
library(tidyverse)
library(monocle3)
library(readxl)
library(here)

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==5) {
  stop("Arguments must be supplied (input file name , metadata, output filename raw, seed, sampling rate).n", call.=FALSE)
} 

# seed
set.seed(args[4])

fraction = as.numeric(args[5])
print("Subsampling fraction")
print(fraction)
print("Project directory")
print(here())

# load data
pca_dmso_raw <- readRDS(args[1])
print("loaded data")

# Setting a PCA cutoff at 25 (hard coded!)
pca_dmso_subset <- pca_dmso_raw %>% 
  dplyr::select(-(PC26:PC50))

# Loading and preparing metadata
metadata = read_excel(args[2])

pca_metadata <- metadata %>% 
  
  mutate(Line = paste0(donor, tumor)) %>% 
  filter(image_CTG == "Imaging") %>%
  dplyr::select(Line, Plate = barcode, screen_ID) %>% 
  #cleaning data to harmonize different naming schemes
  mutate(Plate = if_else(Line == "D020T01" & screen_ID %in% c("HC1092-9", "HC1092-10"), str_replace(Plate, "D020T01", "D020T02"), Plate)) %>%
  mutate(Line = if_else(Line == "D020T01" & screen_ID %in% c("HC1092-9", "HC1092-10"), "D020T02", Line)) %>%
  #accounting for re-imaging of plates
  separate(Plate, c("start", "end"), sep = 8, remove = FALSE) %>% 
  mutate(end = str_sub(end, 2,-1L)) %>% 
  dplyr::select(-Plate, line = Line, start, end, screen_ID) %>%
  left_join(pca_dmso_subset %>% 
              separate(plate, c("start", "end"), sep = 8, remove = FALSE) %>% 
              mutate(end = str_sub(end, 2,-1L)) , .) %>% 
  dplyr::select(everything(), screen_id = screen_ID, -start, -end) %>%
  # reducing dataset size by an order of magnitude
  sample_frac(fraction)

# Import into Monocle 3
## generating metadata objects
pca_anno_df <- pca_metadata %>% dplyr::select(-(PC1:PC25)) %>% 
  mutate(uuid = paste(plate, well, field, object_id, sep = "_")) %>% 
  mutate(uuid2 = uuid) %>% 
  mutate(size = as.numeric(size)) %>%
  mutate(size_log = log(size)) %>%
  as.data.frame() %>% 
  column_to_rownames("uuid2")
print("prepared metadata")


## combining annotation data with harmonized PCA information
pca_matrix <- pca_metadata %>% dplyr::select((PC1:PC25)) %>% as.data.frame() %>% magrittr::set_rownames(pca_anno_df$uuid) %>% magrittr::set_colnames(c(paste0("PC", c(1:25)))) %>% as.matrix()

cce <- new_cell_data_set(pca_matrix %>% t(),
                         cell_metadata = pca_anno_df)


## I manually inject the PCA compression of the data into the object
reducedDims(cce)$PCA <- pca_matrix

print("running gridsearch")
for(i_dist in c(0.01, 0.05, 0.1, 0.5))
  {
  for(j_nn in c(5, 15, 30, 50, 100))
    {
    for(k_res in c(1e-8, 1e-7, 1e-6, 1e-5)){
      print(i_dist)
      print(j_nn)
      print(k_res)
      ## Run UMAP embedding
      cce <- reduce_dimension(cce,
                              reduction_method = "UMAP",
                              umap.min_dist = i_dist,
                              umap.n_neighbors = j_nn,
                              umap.fast_sgd=TRUE, 
                              cores = parallel::detectCores(),
                              verbose = TRUE)
      cce = cluster_cells(cce, resolution=k_res)
      # Save overview plot
      plot_cells(cce) + ggsave(here(paste0("reports/figures/organoid_umap", "_", i_dist, "_", j_nn, "_", k_res,".pdf")))
      # Save result
      saveRDS(cce, paste0(args[3], "_", i_dist, "_", j_nn,  "_", k_res,".Rds"))
      }
  }
}