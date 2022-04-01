# defining lib path
args = commandArgs(trailingOnly=TRUE)
if (args[1] == "remote") {
  .libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/4.0")
  print(.libPaths())
} 

## packages
library(SummarizedExperiment)
library(limma)
library(reshape2)
library(tidyverse)
library(here)
library(DESeq2)

## input
load(here('data/processed/expression/promise_expr.rda'))

organoid_morphology <- read_delim(here::here("references/imaging/visual_classification_organoids.csv"), ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  dplyr::select(line = organoid, morphology = visual_inspection_v2) %>%
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(morphology = if_else(is.na(morphology), "other", morphology))

organoid_size_fit <- readRDS(here::here("data/processed/morphology/organoid_size.Rds")) %>% 
  filter(!line %in% c('D055T01', 'D020T02', 'D021T01')) %>% 
  #filter(!line %in% c('D055T01','D020T02')) %>% 
  mutate(line = as.character(line)) %>% 
  dplyr::select(line, size = x, rep = replicate) %>% 
  distinct() %>% arrange(line) %>%
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(rep = paste0("r", rep))

print("input loaded")
## annotate phenotype group
solid <- organoid_morphology  %>% 
  filter(morphology == "solid") %>%
  .$line
cystic <- organoid_morphology %>% 
  filter(morphology == "cystic") %>%
  .$line

## adding metadata and filtering the promise_expr object
new_col <- colData(promise_expr) %>% as.data.frame() %>%
  left_join(., left_join(
             organoid_size_fit, 
             organoid_morphology)) %>% 
  dplyr::select(-(line:chip_name))
colData(promise_expr) <- cbind(colData(promise_expr), new_col)

promise_expr_filtered <- promise_expr[,promise_expr$line %in% (organoid_size_fit$line %>% unique())]

## long data frame
promise_long <- assays(promise_expr)$expr %>% 
  as_tibble(rownames = 'probe') %>% 
  pivot_longer(values_to = 'expr', names_to = 'id', -probe) %>%
  left_join(as_tibble(rowData(promise_expr), rownames = 'probe')) %>%
  inner_join(as_tibble(colData(promise_expr), rownames = 'id')) %>%
  dplyr::select(-chip_name)

## adding phenotype information
promise_long <- promise_long %>% 
  mutate(phenotype = ifelse(line %in% solid, 'solid', 
                     ifelse(line %in% cystic, 'cystic', 'other'))) #%>% 
 # filter(phenotype != 'other')

## exclude outlier
promise_long_filtered <- promise_long %>% filter(!line %in% c('D054', 'D055', 'D021'))
#promise_long <- promise_long %>% filter(!line %in% c('D054', 'D055'))

## filtering the expression set for top 10% regulated genes 
means <- rowMeans(assays(promise_expr_filtered)$expr)
sd <- apply(assays(promise_expr_filtered)$expr,1,sd)
cv <- sd/means

top_n = cv %>% length()
top_n = top_n * 0.10
keep_probes = cv %>% sort(decreasing = TRUE) %>% .[1:top_n] 
thresh_n = min(keep_probes)

print("probes retained in top variance table")
print(top_n)

promise_expr_filtered_top = promise_expr_filtered[rowData(promise_expr_filtered) %>% rownames() %in% names(keep_probes),]
promise_long_filtered_top = promise_long %>% filter(probe %in% names(keep_probes))

## output
promise_expr_filtered %>% saveRDS(here::here('data/processed/expression/promise_expr_filtered.rds'))
promise_expr_filtered_top %>% saveRDS(here::here('data/processed/expression/promise_expr_filtered_top.rds'))

promise_long %>% saveRDS(here::here('data/processed/expression/promise_expr_tidy.rds'))
promise_long_filtered %>% saveRDS(here::here('data/processed/expression/promise_expr_filtered_tidy.rds'))
promise_long_filtered_top %>% saveRDS(here::here('data/processed/expression/promise_expr_filtered_tidy_top.rds'))

