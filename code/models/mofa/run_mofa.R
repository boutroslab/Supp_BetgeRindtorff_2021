library(here)
library(tidyverse)
library(MOFA2)
library(data.table)


mofa_df <- read_rds(here("data/processed/PhenotypeSpectrum/mofa_df.Rds"))

# brute force reduction
set.seed(123)
mofa_df_filtered <- mofa_df %>% 
  mutate(sample = paste0(sample, "_", group)) %>% 
  mutate(feature = paste0(feature, "_", view)) %>% 
  # log transforming size
  mutate(value = ifelse(feature == "x.0.s.area_expected", log(value), value)) %>%
  filter(view != "mutation") %>%
  # scaling features manually
  group_by(feature, view) %>% 
  mutate(value = if_else(view == "morphology", scale(value), scale(value, center = FALSE))) %>%
  mutate(value = scale(value)) %>%
  ungroup() 

mofa_df_filtered %>% write_csv(here("data/processed/PhenotypeSpectrum/mofa_df_filtered_active_drugs_unscaled.csv"))

input_df <- mofa_df_filtered %>%
  #filter(grepl(sample, pattern = "DMSO")) %>%
  # mutate(feature = paste0(feature, "_", view)) %>% 
  dplyr::select(-group) %>% 
  as.data.table() %>% 
  drop_na()

# printing a diagnostic
input_df %>% sample_frac(0.1) %>% ggplot(aes(value)) + geom_density() + cowplot::theme_cowplot() + facet_wrap(~ view) + ggsave(here("mofa_diagnostic.pdf"))

print("wrangled data")
MOFAobject <- create_mofa(input_df, verbose = TRUE)
print("created MOFA object")
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views = FALSE # default is TRUE
print(data_opts)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors = 5 # default is 15
print(model_opts)

train_opts <- get_default_training_options(MOFAobject)
train_opts$verbose = TRUE
train_opts$maxiter = 10000 # 1000 is default
train_opts$stochastic <- FALSE # default FALSE
train_opts$save_interrupted <- TRUE
print(train_opts)

stochastic_opts <- get_default_stochastic_options(MOFAobject)
print(stochastic_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  #stochastic_options = stochastic_opts, # comment out if not running stochastic inference
  training_options = train_opts
)


outfile = file.path(here("data/processed/PhenotypeSpectrum/mofa_model_active_drugs_unscaled.hdf5"))
# "unscaled" is the result of MOFA without default scaling during the data processing. That said, the data is manually scaled.
# "dmso" is the result of MOFA when run only with DMSO treated organoid data
# "test" contains both "dmso" and "trametinib" treated organoid data
MOFAobject.trained <- run_mofa(MOFAobject, outfile)

