# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/4.0")
print(.libPaths())

# libraries
library(tidyr)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(here)
#library(ggrastr)
#library(cowplot)

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==3) {
  stop("Arguments must be supplied (drug name, fraction_sampling, seed).n", call.=FALSE)
} 
# e.g. Paclitaxel

print(args)


# set seed
set.seed(args[3])

# read complete data
umap_tidy <- read_rds(here::here("data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_tidy.Rds"))

# run queries on data
umap_tidy_selected <- umap_tidy %>% filter(drug == "DMSO" | drug %in% args[1]) %>%
    sample_frac(size = as.numeric(args[2]),
              replace = FALSE)

# save output
filename = str_replace_all(args[1], pattern = " ", "") %>% str_replace_all(., pattern = "/", "_")
umap_tidy_selected %>% write_rds(here::here(paste0("data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_tidy_", filename, ".Rds")))


