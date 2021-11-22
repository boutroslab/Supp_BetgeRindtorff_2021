#!/usr/bin/Rscript
# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

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

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==5) {
  stop("Arguments must be supplied (input file name , output filename harmony, output filename raw, seed, metadata).n", call.=FALSE)
} 

# seed
set.seed(args[4])

# load data
pca_dmso_raw <- readRDS(args[1])

# Setting a PCA cutoff at 25 (hard coded!)
pca_dmso_subset <- pca_dmso_raw %>% 
  dplyr::select(-(PC26:PC50))

# Removing batch effects with Harmony
metadata = read_excel(args[5])

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
  dplyr::select(everything(), screen_id = screen_ID, -start, -end)

# Running Harmony on batches
pca_harmony <- pca_metadata %>% dplyr::select(contains("PC"))
metadata_harmony <- pca_metadata %>% dplyr::select(screen_id, line)

harmony_id <- HarmonyMatrix(
  pca_harmony, metadata_harmony, c("screen_id"),
  do_pca = FALSE,
  verbose = TRUE,
  nclust = 5, # dead, cystic, solid, intermediate, other
  return_object = TRUE
)

# Import into Monocle 3
## generating metadata objects
pca_anno_df <- pca_metadata %>% dplyr::select(-(PC1:PC25)) %>% 
  mutate(uuid = paste(plate, well, field, object_id, sep = "_")) %>% 
  mutate(uuid2 = uuid) %>% 
  mutate(size = as.numeric(size)) %>%
  mutate(size_log = log(size)) %>%
  as.data.frame() %>% 
  column_to_rownames("uuid2")



## combining annotation data with harmonized PCA information
pca_matrix <- pca_metadata %>% dplyr::select((PC1:PC25)) %>% as.data.frame() %>% magrittr::set_rownames(pca_anno_df$uuid) %>% magrittr::set_colnames(c(paste0("PC", c(1:25)))) %>% as.matrix()
pca_matrix_harmony <- harmony_id$Z_corr %>% t() %>% magrittr::set_rownames(pca_anno_df$uuid) %>% magrittr::set_colnames(c(paste0("PC", c(1:25)))) %>% as.matrix()

ods <- new_cell_data_set(pca_matrix_harmony %>% t(),
                         cell_metadata = pca_anno_df)
cce <- new_cell_data_set(pca_matrix %>% t(),
                         cell_metadata = pca_anno_df)

print("included lines:")
print(pca_metadata$line %>% unique())

## I manually inject the PCA compression of the data into the object
reducedDims(ods)$PCA <- pca_matrix_harmony
reducedDims(cce)$PCA <- pca_matrix

# save intermediate result
saveRDS(ods, args[2])
saveRDS(cce, args[3])
print("saved monocle objects at:")
print(args[2])
print(args[3])

