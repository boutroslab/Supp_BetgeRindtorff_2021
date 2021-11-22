
# packages
library(tidyverse)
library(MOFA2)
library(here)
#library(MOFAdata)
library(biomaRt)

# load input 
model <- load_model(here("models/mofa/model.hdf5"))

umap_df <- readRDS(here::here("data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_sampled.Rds"))
promise_long_filtered_top <- readRDS(here::here('data/processed/expression/promise_expr_filtered_tidy_top.rds'))

weights <- model@expectations$Z$single_group %>% as.data.frame() %>% rownames_to_column("id") %>% as_tibble() %>% janitor::clean_names()

drug_anno <- readxl::read_excel(here::here('references/layouts/Compound_Annotation_Libraries_New.xlsx')) %>% distinct(drug = `Compound name`, target = `Primary Target`)

expr_long_gsk <- readRDS(here::here("data/processed/expression/gsk_expr_tidy.rds"))
expr_long_mek <- readRDS(here::here("data/processed/expression/mek_expr_tidy.rds"))

new_morphology <- readRDS(here::here("data/processed/PhenotypeSpectrum/pca_absolute_all_drugs_sampled.Rds")) %>% 
  # dplyr::filter(drug %in% c("DMSO") | (drug %in% c("WYE-125132 (WYE-132)",
  #                                                  "Trametinib (GSK1120212)") & 
  #                                       line %in% c("D046T01", "D018T01"))) %>% 
  # dplyr::filter(drug %in% c("DMSO") | drug %in% c("WYE-125132 (WYE-132)",
  #                                                  "Trametinib (GSK1120212)",
  #                                                 "AZD3759")) %>%
  #dplyr::filter(drug %in% c("DMSO")) %>%
  filter(concentration == "nan") %>%
  mutate(rep = paste0("r", replicate)) %>% 
  mutate(line = substr(line, 1, 4)) %>% 
  mutate(line = paste0(line, "_", rep)) %>%
  # I treat every drug induced phenotype as a dedicated line
  mutate(line = if_else(drug == "DMSO", line, paste0(line, "_", drug))) %>%
  #filter(size_log > 7.5) %>%
  group_by(line) %>% 
  summarise_at(vars(contains("pc")), funs(mean)) %>% 
  ungroup() %>% column_to_rownames("line") %>% as.matrix() %>% 
  #scale(center = TRUE, scale = FALSE) %>% 
  #sweep(2, model@intercepts$morphology_view %>% unlist() %>% unname()) %>% 
  scale(scale = TRUE)

# mapping morphology -> factor 
W <- get_weights(model, "morphology_view")[[1]]
Winv <- pracma::pinv(W)
z_projected <-  new_morphology %*% t(Winv) %>% as.data.frame()%>% rownames_to_column("id")
z_projected_anno <- z_projected %>% 
  mutate(target = substr(id, 9, nchar(id))) %>% 
  mutate(line = substr(id, 1, 7)) %>%
  mutate(target = ifelse(target == "", "DMSO", target))
z_projected_anno_tidy <- z_projected_anno %>% 
  rename(drug = target) %>%
  left_join(drug_anno)
z_projected_anno_test <- z_projected_anno_tidy %>% semi_join(z_projected_anno_tidy %>% 
                                                               dplyr::select(target, drug) %>%
                                                               distinct() %>%
                                            dplyr::count(target) %>% 
                                            filter(n >=5)) %>% 
  filter(drug != "DMSO")

# testing significant movements
drug_projection_test_result <- rbind(
  lm(V1 ~ target, data = z_projected_anno_test) %>% summary() %>% broom::tidy() %>% mutate(factor = "factor1"),
  lm(V2 ~ target, data = z_projected_anno_test) %>% summary() %>% broom::tidy() %>% mutate(factor = "factor2"),
  lm(V3 ~ target, data = z_projected_anno_test) %>% summary() %>% broom::tidy() %>% mutate(factor = "factor3")) 

drug_projection_test_result <- drug_projection_test_result %>% 
  filter(term != "(Intercept)") %>% 
  dplyr::select(term, statistic, factor, p.value)

# saving output
z_projected_anno_tidy %>% saveRDS(here::here("data/morphology/drug_inference.Rds"))
drug_projection_test_result %>% saveRDS(here::here("data/morphology/drug_test.Rds"))
