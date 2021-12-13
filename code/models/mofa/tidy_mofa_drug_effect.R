# packages
library(MOFA2)
library(tidyverse)
#library(here)

set.seed(123)

# input 
## organoid gene expression
promise_long_filtered_top <- readRDS(here::here('data/processed/expression/promise_expr_filtered_tidy_top.rds'))

## organoid size
organoid_size_fit <- readRDS(here::here("data/processed/morphology/organoid_size.Rds")) %>% 
  filter(!line %in% c('D055T01', 'D020T02', 'D021T01')) %>% 
  #filter(!line %in% c('D055T01','D020T02')) %>% 
  mutate(line = as.character(line)) %>% 
  dplyr::select(line, size = x, rep = replicate) %>% 
  distinct() %>% arrange(line) %>%
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(rep = paste0("r", rep))

## organoid morphology
mofa_morphology <- readRDS(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_sampled.Rds")) %>% 
  dplyr::filter(drug %in% c("DMSO") | (drug %in% c("WYE-125132 (WYE-132)") & line == "D046T01")) %>% 
  mutate(rep = paste0("r", replicate)) %>% 
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(line = paste0(line, "_", rep)) %>%
  # I treat every drug induced phenotype as a dedicated line
  mutate(line = if_else(drug == "DMSO", line, paste0(line, "_", drug))) %>%
  #filter(size_log > 7.5) %>%
  group_by(line) %>% 
  summarise_at(vars(contains("pc")), funs(mean)) %>% 
  ungroup() %>% 
  gather(pca, value, -line) %>% rename(feature = pca, sample = line) %>% mutate(view = "morphology_view") %>% 
  mutate(feature = paste0(feature, "_", view))

## organoid drug activity
data('aucroc', package = 'SCOPEAnalysis')
mofa_activity <- aucroc %>% filter(!line %in% c('D055T01', 'D021T01', 'D054T01')) %>%
    mutate(line = substr(line, 1, 4)) %>% 
    expand_grid(., rep = c("r1", "r2")) %>% 
    mutate(sample = paste0(line, "_", rep), value = auroc, feature = paste0("auroc_", drug)) %>% 
    mutate(view = "drug_activity") %>% dplyr::select(sample, feature, value, view) %>% 
    # I average one more time over drug activity per line. This is necessary as D020 was imaged twice. In this particular case, I am creating the average activity score
    group_by(sample, feature, view) %>% summarise(value = mean(value))

# organoid mutation data
mofa_genetics <- read_delim(here::here("data/processed/mutation/Table-S2_Mutations_PDOs_RevisionV4.csv"), delim = ";") %>% 
    janitor::clean_names() %>% 
    mutate(sample = substr(sample, 2,nchar(sample)-1)) %>% 
    mutate(sample = paste0("D", sample)) %>% 
    dplyr::filter(!sample %in% c("D021", "D015")) %>%
    expand_grid(replicate = c("r1", "r2")) %>% 
    mutate(sample = paste0(sample, "_", replicate)) %>% 
    dplyr::select(sample, feature = symbol, everything()) %>% dplyr::select(sample, feature) %>% 
    mutate(value = 1) %>%
    complete(sample, feature, fill = list(value = 0)) %>% 
    distinct(sample, feature, value) %>%
    mutate(view = "mutation")

# define MOFA object and run MOFA
mofa_size <- organoid_size_fit %>% 
  mutate(line = paste0(line, "_", rep)) %>%
  group_by(line) %>% summarise(value = mean(size)) %>% 
  mutate(feature = "size",
         view = "size_view") %>% 
  drop_na() %>% 
  dplyr::select(sample = line, feature, view, value)

mofa_expression <- promise_long_filtered_top %>% 
  # renaming columns
  mutate(line = paste0(line, "_", rep)) %>%
  dplyr::select(sample = line,
                feature = symbol, # setting feature to symbol
                value = expr) %>% 
  mutate(view = "expression") %>% 
  # averaging feature value
  group_by(feature, sample, view) %>%
  summarise(value = mean(value)) %>%
  # renaming feature jic
  mutate(feature = paste0(feature, "_", view)) %>%
  dplyr::select(sample, feature, view, value) %>% 
  drop_na() %>% 
  filter(feature != "_expression")

mofa_classification <- promise_long_filtered_top %>% 
  dplyr::select(sample = line, value = morphology) %>% 
  mutate(view = "view_classification", feature = "classification") %>% 
  mutate(value = factor(value) %>% as.numeric(),
         value = value -1,
         value = as.logical(value)) %>%
  distinct()

input_df = rbind(mofa_size,
                 mofa_morphology,
                 #mofa_classification,
                 mofa_activity,
                 mofa_genetics,
                 mofa_expression
                 ) %>% 
  data.table::as.data.table()

MOFAobject <- create_mofa(input_df, verbose = TRUE)
print("created MOFA object")
plot_data_overview(MOFAobject) #%>% ggsave(here::here("reports/figures/mofa_object.pdf"))
system("mv Rplots.pdf reports/figures/mofa_object.pdf")

## setting option
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views = TRUE # default is TRUE
print(data_opts)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors = 3 # default is 15
print(model_opts)

train_opts <- get_default_training_options(MOFAobject)
train_opts$verbose = TRUE
train_opts$maxiter = 1000 # 1000 is default
train_opts$stochastic <- FALSE # default FALSE
train_opts$save_interrupted <- TRUE
print(train_opts)

stochastic_opts <- get_default_stochastic_options(MOFAobject)
print(stochastic_opts)

# running model
outfile = file.path(here::here("models/mofa/model_drug_effect.hdf5"))
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  #stochastic_options = stochastic_opts, # comment out if not running stochastic inference
  training_options = train_opts
)
MOFAobject.trained <- run_mofa(MOFAobject, outfile)