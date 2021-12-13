# example
# Rscript --vanilla src/models/umap/annotate_size.R data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_sampled.Rds FALSE
# Rscript --vanilla src/models/umap/annotate_size.R data/processed/PhenotypeSpectrum/umap_absolute_all_drugs.Rds TRUE

# input
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==2) {
  stop("Arguments must be supplied (input file name, remote? boolean).n", call.=FALSE)
} 

print(args)

if (args[2] == TRUE) {
  # defining lib path
  .libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
  print("running in remote mode, make sure you are running code in R/3.6.0")
  print(.libPaths())
}

## library 
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(readr)
library(here)
library(fitdistrplus)

## loading data
print("loading data")
umap_df <- read_rds(here::here(args[1]))

# ## overall size
# line_param <- umap_df %>% #filter(partition %in% c(1,2)) %>% 
#   nest(-line, -replicate) %>% 
#   mutate(fit = map(data, ~ fitdistrplus::fitdist(.x$size, "lnorm")),
#          param = map(fit, ~ .x$estimate %>% broom::tidy()))

# df <- line_param %>% unnest(param) %>% 
#   filter(names == "meanlog") %>% 
#   group_by(line) %>% 
#   mutate(mean_meanlog = mean(x)) %>% 
#   arrange(mean_meanlog) %>% 
#   ungroup() %>%
#   mutate(line = factor(line, levels = .$line %>% unique()))

# organoid_size_fit <- df %>% dplyr::select(line, replicate, names, x, mean_meanlog)
# organoid_size_fit %>% saveRDS(here::here("data/processed/morphology/organoid_size_fit.Rds"))

# ### no fit
# df <- umap_df %>% 
#   group_by(line, replicate) %>% 
#   summarise(mean = mean(size_log),
#             median = median(size_log))

# df %>% saveRDS(here::here("data/processed/morphology/organoid_size.Rds"))


# ## DMSO size
# line_param <- umap_df %>% #filter(partition %in% c(1,2)) %>% 
#   filter(drug == "DMSO") %>%
#   nest(-line, -replicate) %>% 
#   mutate(fit = map(data, ~ fitdistrplus::fitdist(.x$size, "lnorm")),
#          param = map(fit, ~ .x$estimate %>% broom::tidy()))

# df <- line_param %>% unnest(param) %>% 
#   filter(names == "meanlog") %>% 
#   group_by(line) %>% 
#   mutate(mean_meanlog = mean(x)) %>% 
#   arrange(mean_meanlog) %>% 
#   ungroup() %>%
#   mutate(line = factor(line, levels = .$line %>% unique()))

# organoid_size_fit <- df %>% dplyr::select(line, replicate, names, x, mean_meanlog)
# organoid_size_fit %>% saveRDS(here::here("data/processed/morphology/organoid_size_DMSO_fit.Rds"))

# ### no fit
# df <- umap_df %>% 
#   filter(drug == "DMSO") %>%
#   group_by(line, replicate) %>% 
#   summarise(mean = mean(size_log),
#             median = median(size_log))

# df %>% saveRDS(here::here("data/processed/morphology/organoid_size_DMSO.Rds"))


# ## by drug size - CAUTION make sure to run the estimation of drug effects on size with the complete dataset, not a subsample
# line_param <- umap_df %>% #filter(partition %in% c(1,2)) %>% 
#   nest(-line, -replicate, -drug, -concentration) %>% 
#   mutate(fit = map(data, ~ fitdistrplus::fitdist(.x$size, "lnorm")),
#          param = map(fit, ~ .x$estimate %>% broom::tidy()))

# df <- line_param %>% unnest(param) %>% 
#   filter(names == "meanlog") %>% 
#   group_by(line) %>% 
#   mutate(mean_meanlog = mean(x)) %>% 
#   arrange(mean_meanlog) %>% 
#   ungroup() %>%
#   mutate(line = factor(line, levels = .$line %>% unique()))

# organoid_size_fit <- df %>% dplyr::select(line, replicate, names, x, mean_meanlog)
# organoid_size_fit %>% saveRDS(here::here("data/processed/morphology/organoid_size_drug_fit.Rds"))

### no fit
df <- umap_df %>% 
  group_by(line, plate, replicate, well, drug, concentration) %>% 
  summarise(mean = mean(size_log),
            median = median(size_log))

df %>% saveRDS(here::here("data/processed/morphology/organoid_size_drug.Rds"))

## channel intensity
df <- umap_df %>% 
  group_by(line, plate, replicate, well, drug, concentration) %>% 
  summarise(permeability = mean(permeability),
            dapi = mean(dapi),
            actin = mean(actin))

df %>% saveRDS(here::here("data/processed/morphology/organoid_intensity_drug.Rds"))
