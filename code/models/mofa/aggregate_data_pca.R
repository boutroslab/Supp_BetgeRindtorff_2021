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

# read complete data
pca_tidy <- read_rds(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_tidy.Rds"))

# run queries on data
pca_tidy_aggregate <- pca_tidy %>% 
  mutate(concentration = if_else(drug == "DMSO", "nan", concentration)) %>%
  dplyr::group_by(line, drug, concentration, replicate) %>% 
  summarise_at(vars(contains("pc")), funs(mean)) %>% 
  ungroup() %>%
  left_join(pca_tidy %>%
    count(line, drug, concentration, replicate))

# save output
pca_tidy_aggregate %>% write_rds(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_aggregate.Rds"))


