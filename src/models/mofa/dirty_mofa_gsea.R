library(tidyverse)
library(MOFA2)
library(here)
library(MOFAdata)
library(ReactomePA)
library(biomaRt)


model <- load_model(here("data/processed/PhenotypeSpectrum/mofa_model_active_drugs_unscaled.hdf5"))
load(here("notebooks/SCOPEAnalysis/data/promise_expr.rda")) 
feature_anno <- read_csv(here("data/interim/FeatureAnalysis/feature_annotation.csv"))

PATH = "/dkfz/groups/shared/OE0049/B110-Isilon2/promise/"

umap_tidy <- read_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_tidy.Rds"))
umap_sampled <- read_rds(paste0(PATH, "data/processed/PhenotypeSpectrum/umap_absolute_all_drugs_sampled.Rds"))

# fomrmatting model output
weights <- model@expectations$Z$single_group %>% as.data.frame() %>% rownames_to_column("id") %>% as_tibble() %>% janitor::clean_names()
loadings_morphology <- model@expectations$W$morphology %>% as.data.frame() %>% rownames_to_column("id") %>% as_tibble() %>% janitor::clean_names()
loadings_expression <- model@expectations$W$expression %>% as.data.frame() %>% rownames_to_column("id") %>% as_tibble() %>% janitor::clean_names()
intercepts_morphology = model@intercepts$morphology$single_group %>% as.data.frame() %>% rownames_to_column("id") %>% as_tibble() %>% janitor::clean_names()
intercepts_expression = model@intercepts$morphology$single_group %>% as.data.frame() %>% rownames_to_column("id") %>% as_tibble() %>% janitor::clean_names()

# converting gene names
set.seed(1343)

model@expectations$W$expression

feature_metadata_conversion <- model@features_metadata %>% 
  mutate(gene = substr(feature, 1, nchar(feature)-11))

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
dict= biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = feature_metadata_conversion$gene, bmHeader = T, mart = mart)
colnames(dict)=c("ensembl_gene_id","hgnc_symbol")

# TODO In an unorthodox way, I arbitrarily sample a ENSG id linking to each hgnc symbol, by taking the first ENSG in cases where there are many
dict_sampled <- dict %>% 
  nest(-hgnc_symbol) %>% 
  mutate(short = purrr::map(data, ~ .x %>% head(1))) %>% 
  dplyr::select(-data) %>% 
  unnest() 

df <- dict_sampled %>% 
  mutate(id = paste0(hgnc_symbol, "_expression")) %>%
  semi_join(loadings_expression)

# TODO check dropping of genes in merge = the reason I introduced a semi_join above
model_ensg <- model
index_model <- rownames(model_ensg@expectations$W$expression) %in% df$id

tmp = model_ensg@expectations$W$expression[index_model,]
rownames(tmp) = df$ensembl_gene_id
model_ensg@expectations$W$expression <- tmp


model_ensg



gene_sets = msigdbr::msigdbr(#category = "C2",
  species = "Homo sapiens")

gene_sets_matrix <- gene_sets %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  mutate(present = 1,
         gene_symbol = paste0(gene_symbol, "_expression")) %>% 
  pivot_wider(names_from = gene_symbol,
              values_from = present,
              values_fill = 0) %>%
  as.data.frame() %>% column_to_rownames("gs_name") %>% as.matrix()

enrichment.parametric_all <- MOFA2::run_enrichment(model,
                                                   view = "expression", factors = 1:5,
                                                   feature.sets = gene_sets_matrix,
                                                   sign = "all",
                                                   statistical.test = "parametric"
)

enrichment.parametric_negative <- MOFA2::run_enrichment(model,
                                                        view = "expression", factors = 1:5,
                                                        feature.sets = gene_sets_matrix,
                                                        sign = "negative",
                                                        statistical.test = "parametric"
)

enrichment.parametric_positive <- MOFA2::run_enrichment(model,
                                                        view = "expression", factors = 1:5,
                                                        feature.sets = gene_sets_matrix,
                                                        sign = "positive",
                                                        statistical.test = "parametric"
)

enrichment.cor.adj.parametric_all <- MOFA2::run_enrichment(model,
                                                           view = "expression", factors = 1:5,
                                                           feature.sets = gene_sets_matrix,
                                                           sign = "all",
                                                           statistical.test = "cor.adj.parametric"
)

enrichment.cor.adj.parametric_negative <- MOFA2::run_enrichment(model,
                                                                view = "expression", factors = 1:5,
                                                                feature.sets = gene_sets_matrix,
                                                                sign = "negative",
                                                                statistical.test = "cor.adj.parametric"
)

enrichment.cor.adj.parametric_positive <- MOFA2::run_enrichment(model,
                                                                view = "expression", factors = 1:5,
                                                                feature.sets = gene_sets_matrix,
                                                                sign = "positive",
                                                                statistical.test = "cor.adj.parametric"
)

enrichment.permutation_all <- MOFA2::run_enrichment(model,
                                                    view = "expression", factors = 1:5,
                                                    feature.sets = gene_sets_matrix,
                                                    sign = "all",
                                                    statistical.test = "permutation"
)

enrichment.permutation_negative <- MOFA2::run_enrichment(model,
                                                         view = "expression", factors = 1:5,
                                                         feature.sets = gene_sets_matrix,
                                                         sign = "negative",
                                                         statistical.test = "permutation"
)

enrichment.permutation_positive <- MOFA2::run_enrichment(model,
                                                         view = "expression", factors = 1:5,
                                                         feature.sets = gene_sets_matrix,
                                                         sign = "positive",
                                                         statistical.test = "permutation"
)

enrichment.parametric_all %>% saveRDS(here("data/processed/PhenotypeSpectrum/mofa_gsea_parametric_all.Rds"))
saveRDS(enrichment.parametric_negative, here("data/processed/PhenotypeSpectrum/mofa_gsea_parametric_down.Rds"))
saveRDS(enrichment.parametric_positive, here("data/processed/PhenotypeSpectrum/mofa_gsea_parametric_up.Rds"))

enrichment.cor.adj.parametric_all %>% saveRDS(here("data/processed/PhenotypeSpectrum/mofa_gsea_adjparametric_all.Rds"))
saveRDS(enrichment.cor.adj.parametric_negative, here("data/processed/PhenotypeSpectrum/mofa_gsea_adjparametric_down.Rds"))
saveRDS(enrichment.cor.adj.parametric_positive, here("data/processed/PhenotypeSpectrum/mofa_gsea_adjparametric_up.Rds"))

enrichment.permutation_all %>% saveRDS(here("data/processed/PhenotypeSpectrum/mofa_gsea_permutation_all.Rds"))
saveRDS(enrichment.permutation_negative, here("data/processed/PhenotypeSpectrum/mofa_gsea_permutation_down.Rds"))
saveRDS(enrichment.permutation_positive, here("data/processed/PhenotypeSpectrum/mofa_gsea_permutation_up.Rds"))