# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==5) {
  stop("Arguments must be supplied (raw input file name, harmony input file name, seed, sampling_rate, source for feature intensity function).n", call.=FALSE)
} 

print(args)

# seed
set.seed(args[3])


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

# I wish I could solve my path problems with the here package. 
#PATH = "/dkfz/groups/shared/OE0049/B110-Isilon2/promise/"
PATH = paste0(here::here(), "/")

## TODO the two scripts need to be un-commented for renewed execution of the code
# hacky way - extracting intensity features for objecy creation
source(here::here("src/models/umap/annotate_ldc.R"))
source(here::here("src/models/umap/annotate_features.R"))

annotate_features(here::here(args[5]))

# Loading data and formatting for easier access, commented non-essential code
obj <- readRDS(paste0(PATH, args[1]))
obj_h <- readRDS(paste0(PATH, args[2]))
# dmso_h <- readRDS(paste0(PATH, "data/interim/PhenotypeSpectrum/harmony_umap_absolute_dmso.Rds"))
# dmso <- readRDS(paste0(PATH, "data/interim/PhenotypeSpectrum/hdf5_umap_absolute_dmso.Rds"))

# Formatting UMAP tables for easy access
umap_tidy <- reducedDims(obj)$UMAP %>% cbind(colData(obj)) %>% as_tibble() %>% janitor::clean_names() 
pca_tidy <- reducedDims(obj)$PCA %>% cbind(colData(obj)) %>% as_tibble() %>% janitor::clean_names()

# Formatting other UMAP projections
umap_tidy_h <- reducedDims(obj_h)$UMAP %>% cbind(colData(obj_h)) %>% as_tibble() %>% janitor::clean_names()
pca_tidy_h <- reducedDims(obj_h)$PCA %>% cbind(colData(obj_h)) %>% as_tibble() %>% janitor::clean_names()
# ## Spiking in LDC classifier data to DMSO-only data
# umap_tidy_dmso_h <- reducedDims(dmso_h)$UMAP %>% cbind(colData(dmso_h)) %>% as_tibble() %>% janitor::clean_names()
# umap_tidy_dmso <- reducedDims(dmso)$UMAP %>% cbind(colData(dmso)) %>% as_tibble() %>% janitor::clean_names()

# adding intensity values to the umap/pca_tidy object
intensity <- read_rds(paste0(PATH, "data/interim/FeatureAnalysis/feature_intensity.rds")) %>%
  magrittr::set_colnames(c("actin", "permeability", "dapi"))
# merging the columns
umap_tidy <- cbind(umap_tidy, intensity)
pca_tidy <- cbind(pca_tidy, intensity)

umap_tidy_h <- cbind(umap_tidy_h, intensity)
pca_tidy_h <- cbind(pca_tidy_h, intensity)

# adding cluster and segment-level data
umap_tidy <- cbind(umap_tidy, cluster = clusters(obj), partition = partitions(obj)) %>% as_tibble()
pca_tidy <- cbind(pca_tidy, cluster = clusters(obj), partition = partitions(obj)) %>% as_tibble()

umap_tidy_h <- cbind(umap_tidy_h, cluster = clusters(obj_h), partition = partitions(obj_h)) %>% as_tibble()
pca_tidy_h <- cbind(pca_tidy_h, cluster = clusters(obj_h), partition = partitions(obj_h)) %>% as_tibble()

# TODO adding LDC values
ldc = readRDS(here::here("data/processed/ldc_viability.rds")) %>%
  as.data.frame() %>% dplyr::select(uuid, prediction, prob_dead, prob_live)

umap_tidy <- left_join(umap_tidy, ldc, by = "uuid")
pca_tidy  <- left_join(pca_tidy, ldc, by = "uuid")
umap_tidy_h <- left_join(umap_tidy_h, ldc, by = "uuid")
pca_tidy_h  <- left_join(pca_tidy_h, ldc, by = "uuid")

# checking structure
print(str(umap_tidy))

# I subsample parts of my data. 
umap_sampled <- umap_tidy %>%
  sample_frac(size = as.numeric(args[4]),
              replace = FALSE)

pca_sampled <- pca_tidy %>%
  sample_frac(size = as.numeric(args[4])/2, # pca object is larger and requires lower sampling rate
              replace = FALSE)

umap_sampled_h <- umap_tidy_h %>%
  sample_frac(size = as.numeric(args[4]),
              replace = FALSE)

pca_sampled_h <- pca_tidy_h %>%
  sample_frac(size = as.numeric(args[4])/2,
              replace = FALSE)

# write to file
umap_tidy %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_tidy.Rds"))
pca_tidy %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_tidy.Rds"))
umap_sampled %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_sampled.Rds"))
pca_sampled %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_sampled.Rds"))

# other output data for comparison
umap_tidy_h %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/harmony_umap_absolute_all_drugs_tidy.Rds"))
pca_tidy_h %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/harmony_pca_absolute_all_drugs_tidy.Rds"))
umap_sampled_h %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/harmony_umap_absolute_all_drugs_sampled.Rds"))
pca_sampled_h %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/harmony_pca_absolute_all_drugs_sampled.Rds"))

# umap_tidy_dmso_h %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/harmony_umap_absolute_dmso_tidy.Rds"))
# umap_tidy_dmso %>% write_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/umap_absolute_dmso_tidy.Rds"))
