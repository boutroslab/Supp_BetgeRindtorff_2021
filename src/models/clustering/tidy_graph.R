library(tidyverse)
library(monocle3)
library(here)


length_in <- c(1e-8, 1e-7, 1e-6, 1e-5)
for(i in length_in){
  path = here(paste0("data/processed/PhenotypeSpectrum/louvain_absolute_all_drugs_", i))
  print(path)
  obj <- readRDS(path)
  

  plot_cells(obj) + ggsave(here(paste0("reports/figures/louvain_absolute_all_drugs_", i, ".pdf")))
  plot_cells(obj, color_cells_by="partition", group_cells_by="partition") + ggsave(here(paste0("reports/figures/louvain_absolute_all_drugs_partition", i, ".pdf")))
}