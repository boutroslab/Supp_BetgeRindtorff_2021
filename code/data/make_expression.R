# title: "Generation of baseline expression file for PROMISE manuscript"
# author: "Benedikt Rauscher, Niklas Rindtorff"
# date: "7/2/2021"
# note: refactored from "br_baseline_expression.Rmd"

# checking input and defining lib path
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==2) {
  stop("arguments must be supplied (remote / local, baseline / all).n", call.=FALSE)
} 

if (args[1] == "remote") {
  .libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/4.0")
  print(.libPaths())
} 

## TODO integrate all treatment analysis
if (args[2]== "all"){
  stop("aggregate normalization with baseline and drug treatments has not yet been implemented. aborting.")
}

# loading libraries
library(tidyverse)
library(affy)
library(org.Hs.eg.db)
library(here)
require("biomaRt")

# Data preparation
# ================
# In order to perform subtype classification a number of annotations is necessary that can help to link microarray probes to Entrez IDs required for subtype classification. The following code generate this mapping and in addition creates a list of reference genes whose expression will be required for the subtype classification.

## id mapping derived from ENSEMBL biomart
id_mapping <- read_tsv(here::here('references/microarrays/hugene_mapping.txt'))
id_mapping <- id_mapping %>% `colnames<-`(c('ensg', 'probe', 'symbol', 'dummy')) %>% 
  dplyr::select(-dummy)
x <- org.Hs.egSYMBOL
mapping_all <- as.list(x[mappedkeys(x)])
mapping_all <- tibble(entrez = names(mapping_all), symbol = mapping_all) %>% 
  `colnames<-`(c('entrez', 'symbol')) %>% unnest(symbol)
id_mapping <- id_mapping %>% inner_join(mapping_all) %>% drop_na() %>% distinct()

## U133 plus2 microarray probe annotation
id_mapping_u133 <- read_tsv(here::here('references/microarrays/affy_mapping2.txt')) %>%
  dplyr::select(`Ensembl Gene ID`, `Affy HG U133-PLUS-2 probeset`) %>%
  `colnames<-`(c('ensg', 'probe')) %>%
  inner_join(id_mapping %>% dplyr::select(-probe))

# Raw data normalization
# ======================

## adapt path to point to folder with cel files. 
fn <- list.files(here::here('data/raw/expression/microarray_manuscript/baseline/CEL_files'), recursive = T,
                 full.names=T, pattern = '*CEL') %>% .[!grepl('HuGene', .)]

sn <- paste('sample', 1:length(fn), sep='_')

if (args[2] == "all") {
  # appending other drug treatments
  fn <- c(fn, list.files(here::here('data/raw/expression/microarray_manuscript/GSK_inhibition/CEL_files'), recursive = T,
                 full.names=T, pattern = '*CEL') %>% .[!grepl('HuGene', .)],
              list.files(here::here('data/raw/expression/microarray_manuscript/MEK_inhibition/CEL_files'), recursive = T,
                 full.names=T, pattern = '*CEL') %>% .[!grepl('HuGene', .)]
                 )

  sn <- paste('sample', 1:length(fn), sep='_')
} 

## read chips
orgas_bl <- ReadAffy(filenames = fn, sampleNames = sn)

## calculate RMA corrected expression values
ex_bl <- expresso(orgas_bl,
                  bgcorrect.method='rma', 
                  normalize.method='quantiles',
                  pmcorrect.method='pmonly',
                  summary.method='avgdiff')

## remember the full matrix (for GEO submission)
full_exp_mat <- ex_bl %>% exprs()

## sample object
samples_nb <- tibble(filenames = fn, sample_names = sn) %>% 
  mutate(chip_name = basename(filenames)) %>% 
  tidyr::extract(filenames, 'line', 
          regex = '[PD]*0*(\\d{1,3}).+$', 
          remove=F) %>%
  group_by(line) %>%
  mutate(rep = paste0('r', 1:n())) %>% 
  ungroup() %>%
  mutate(line = ifelse(nchar(line) == 2, paste0('D0', line), paste0('D00', line)),
         id = paste(line, rep, sep='_'))

## matrix
log_expr <- full_exp_mat %>% log()
stopifnot(identical(samples_nb$sample_names, colnames(log_expr)))
colnames(log_expr) <- samples_nb$id

# Generate improved probe annotation using biomart
# ====================================

## pulling probes from biomart

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hg_u133_plus_2",
    "ensembl_gene_id",
    "gene_biotype",
    "entrezgene_id",
    "external_gene_name"),
  filter = "affy_hg_u133_plus_2",
  values = rownames(log_expr), uniqueRows=TRUE)

# formatting
id_mapping_u133_mart <- annotLookup %>% as_tibble() %>% 
  dplyr::select(probe = affy_hg_u133_plus_2, 
                ensg = ensembl_gene_id, 
                entrez = entrezgene_id, 
                symbol = external_gene_name) %>% 
  arrange(probe, entrez, symbol, ensg) %>% distinct() %>% group_by(probe) %>% dplyr::slice(1) %>% ungroup() 

row_data_mart = tibble(probe = rownames(log_expr)) %>% 
  left_join(id_mapping_u133_mart)
stopifnot(identical(row_data_mart$probe, rownames(log_expr)))

print("missing annotation statistic")
print(row_data_mart %>% mutate(missing = is.na(symbol)) %>% dplyr::count(missing))


# Generate summarized experiment object.
# ====================================

## row data
row_data <- row_data_mart
#### manual annotation
# row_data <- tibble(probe = rownames(log_expr)) %>% 
#   left_join(distinct(id_mapping_u133) %>%
#               group_by(probe) %>% dplyr::slice(1) %>% ungroup())
# identical(row_data$probe, rownames(log_expr))

# print("missing annotation statistic")
# print(row_data %>% mutate(missing = is.na(symbol)) %>% dplyr::count(missing))

## col data
col_data <- samples_nb %>% dplyr::select(id, line, rep, chip_name) %>% 
  as.data.frame() %>% column_to_rownames('id')
identical(rownames(col_data), colnames(log_expr))

## combine into summarized experiment object
promise_expr <- SummarizedExperiment::SummarizedExperiment(
  assays = list(expr = log_expr),
  rowData = row_data %>% as.data.frame() %>% column_to_rownames('probe'),
  colData = col_data
)


if (args[2] != "all") {
  save(promise_expr, file = here::here('data/processed/expression/promise_expr.rda'))
} 

if (args[2] == "all") {
  save(promise_expr, file = here::here('data/processed/expression/promise_expr_all.rda'))
} 

# Session info
sessionInfo()


