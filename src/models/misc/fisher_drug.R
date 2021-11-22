library(tidyverse)
library(here)

df_cluster <- readRDS(here("data/clustered_umap.Rds"))

chisq_drug <- df_cluster %>% 
  dplyr::select(mclust_complete_classification, drug) %>%
  filter(drug != "DMSO") %>%
  filter(drug != "Staurosporine_500nM") %>%
  mutate(mclust_complete_classification = as.numeric(mclust_complete_classification)) %>%
  group_by(mclust_complete_classification, drug) %>% 
  summarize(n = n()) %>% 
  spread(drug, n, fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("mclust_complete_classification") %>% 
  as.matrix() %>%
  #chisq.test()
  stats::fisher.test(simulate.p.value = TRUE)

saveRDS(chisq_drug, here("data/drug_contig/fisher.Rds"))