# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==5) {
  stop("Arguments must be supplied.n", call.=FALSE)
} 

print(args)


# libraries
library(tidyr)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(rhdf5)
library(here)

pca_path <- here::here(args[1])
h5ls(pca_path)
hdf5_pca <- h5read(pca_path, "features_organoids") %>% as_tibble()
hdf5_pca_meta <- h5read(pca_path, "metadata_organoids") %>% as.matrix %>% t() %>% as_tibble() %>% magrittr::set_colnames(h5read(pca_path, "metadata_names_organoids"))

print(hdf5_pca_meta$Line %>% table())

hdf5_pca <- cbind(hdf5_pca_meta, hdf5_pca) %>% 
  arrange(Line, Plate, Well) %>%
  janitor::clean_names() %>%
  magrittr::set_colnames(colnames(.) %>% str_replace('v', 'PC'))
  
print(head(hdf5_pca))

saveRDS(hdf5_pca, here(args[2]))

print(here(args[2]))

## processing PCA metadata
variance_explained <- read_csv(here::here(args[3])) %>% 
  dplyr::select(PC = X1, eigenvalue = `0`) %>% mutate(var_explained = eigenvalue/sum(eigenvalue)) %>% 
  mutate(sum_explained = cumsum(var_explained)) %>% mutate(PC = PC +1)

components <- read_csv(here::here(args[4]))
feature_annotation <- read_csv(here::here(args[5])) 
annotated_PCA <- components %>% dplyr::select(-X1) %>% magrittr::set_colnames(paste0("PC", 1:ncol(.))) %>% cbind(feature_annotation, .) %>% as_tibble()

# writing PCA annotation to disk
variance_explained %>% saveRDS(here::here("data/processed/PhenotypeSpectrum/pca_variance.Rds"))
annotated_PCA %>% saveRDS(here::here("data/processed/PhenotypeSpectrum/pca_loading.Rds"))



