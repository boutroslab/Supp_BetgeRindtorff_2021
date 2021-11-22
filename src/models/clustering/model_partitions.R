# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

library(nnet)
library(tidyverse)
library(tidyr)
library(here)

umap_df

# preparing data matrix and using partition 1 as reference level
df <- umap_df %>% 
  dplyr::select(drug, line, partition, plate, screen_id) %>% 
  mutate(drug = factor(drug),
         line = factor(line),
         plate = factor(plate),
         screen_id = factor(screen_id))

df$partition %>% table()

model_intercept = multinom(partition ~ 1, data = df, model = TRUE)
model_line = multinom(partition ~ line, data = df, model = TRUE)
model_plate = multinom(partition ~ plate, data = df, model = TRUE)
model_screenid = multinom(partition ~ screen_id, data = df, model = TRUE)
#model_drug = multinom(partition ~ drug, data = df, model = TRUE)

aic_multinomial <- AIC(model_intercept, model_line, model_plate, model_screenid) %>% rownames_to_column("model") %>% arrange(AIC) 

# writing
aic_multinomial %>% write_csv(here::here("reports/tables/model_partition.csv"))
model_intercept %>% write_rds(here::here("models/cluster/model_intercept.rds"))
model_line %>% write_rds(here::here("models/cluster/model_line.rds"))
model_plate  %>% write_rds(here::here("models/cluster/model_plate.rds"))
model_screenid %>% write_rds(here::here("models/cluster/model_screenid.rds"))
