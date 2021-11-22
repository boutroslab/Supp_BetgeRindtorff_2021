#!/usr/bin/Rscript

# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

# packages
library(tidyr)
library(readr)
library(dplyr)
library(magrittr)
library(monocle3)
library(here)

## TODO reenable command line knobs after refactoring and removing script from tidy_umap.R
# run:
#args = commandArgs(trailingOnly=TRUE)

# not run: 
# args <- c("data/interim/FeatureAnalysis/organoid_viability/human/results", "data/interim/PhenotypeSpectrum/harmony_umap_absolute_all_drugs.Rds")

#print(args)

#if (!length(args)==2) {
#  stop("Arguments must be supplied (input directory logs, input monocle object).n", call.=FALSE)
#} 

# hacky way to define input
ldc_input = "data/interim/FeatureAnalysis/organoid_viability/human/results" # args[1]
monocle_input = "data/interim/PhenotypeSpectrum/archive_pre_qc/hdf5_umap_absolute_all_drugs.Rds" # args[2]

#drop_lines = c("D054T01", "D055T01", "D046T01", "D020T02", "D010T01")
#drop_plates = c("D027T01P906L03", "D020T01P906L03", "D013T01P001L02")

# find classification_logs and aggregate results for UMAP filtering
log <- list.files(ldc_input, full.names = TRUE, pattern = "log.csv") %>% # the log files contain the granular classification info
  as_tibble() %>% janitor::clean_names() %>%
  # reading logs
  mutate(data = purrr::map(value, ~ .x %>% read_csv)) %>% 
  # crunching filenames to extract organoid line name in a hacky way
  separate(value, into = c("path", "ending"), sep = "_", fill = "right") %>% 
  mutate(line = substr(ending, nchar(ending)-6, nchar(ending))) %>%
  dplyr::select(everything(), -path, -ending) %>% 
  # counting the number of objects
  mutate(n_log = purrr::map(data, ~ .x %>% nrow())) %>% 
  unnest(n_log)

colnames(log)


# accessing the monocle object
obj <- readRDS(monocle_input)
umap_tidy <- reducedDims(obj)$UMAP %>% cbind(colData(obj)) %>% as_tibble() %>% janitor::clean_names()
  

# TEST: the number of objects in the log file has to match the number of objects in the monocle object.
obj_n <- umap_tidy %>% dplyr::count(line)
ldc_n <- log %>% dplyr::select(line, n = n_log)
stopifnot(obj_n == ldc_n)

# I am saving my final result
## monocle object
# write object with added LDC information.
df <- log %>% dplyr::select(everything(), -n_log) %>% 
  unnest(data) %>% 
  dplyr::select(-line) %>% janitor::clean_names() %>% 
  cbind(colData(obj), .) 
#pData(obj) <- df
#obj %>% saveRDS(monocle_input)

## dedicated ldc object
df %>% write_rds(here::here("data/processed/ldc_viability.rds"))