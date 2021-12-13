library(here)
library(tidyverse)
library(readxl)
library(SummarizedExperiment)

# load assets
print("loading gene expression data")
load(here("data/processed/expression/promise_expr.rda"))

print("loading morphology annotation")
organoid_morphology <- read_delim(here::here("references/imaging/visual_classification_organoids.csv"), ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  dplyr::select(line = organoid, morphology = visual_inspection_v2)

print("loading imaging metadata")
morph_metadata = read_excel("references/Screenings_Imaging_Manuscript.xlsx")

## annotate phenotype group
cystic_l <- organoid_morphology %>% filter(morphology == "cystic") %>%.$line 
dense_l <- organoid_morphology %>% filter(morphology == "solid") %>%.$line 

## TODO simplify



# gene expression
## long data frame
promise_long <- assays(promise_expr)$expr %>% 
  as_tibble(rownames = 'probe') %>% 
  pivot_longer(values_to = 'expr', names_to = 'id', -probe) %>%
  left_join(as_tibble(rowData(promise_expr), rownames = 'probe')) %>%
  inner_join(as_tibble(colData(promise_expr), rownames = 'id')) %>%
  select(-chip_name)

## adding phenotype information
promise_long <- promise_long %>% 
  mutate(phenotype = ifelse(line %in% dense_l, 'solid', 
                      ifelse(line %in% cystic_l, 'cystic', 'other')))

## exclude outlier - using a leftjoin for that
# promise_long <- promise_long %>% filter(!line %in% c('D054', 'D055', 'D021')) %>% 
#   filter(!line %in% c('D046', 'D010'))

## select most highly expressed probe to represent each gene
select_probes <- promise_long %>% group_by(symbol, probe) %>% 
  summarise(avg_probe = mean(expr)) %>% ungroup() %>%
  group_by(symbol) %>% top_n(1, avg_probe) %>% ungroup() %>% pull(probe)

## summarize replicates
gene_expr <- promise_long %>% 
  # group_by(line, symbol, probe, phenotype) %>%
  # summarise(expr = mean(expr)) %>% ungroup() %>%
  filter(probe %in% select_probes)

# imaging features
print('loading imaging data')

# TODO #118 integrate Common_Feautres_human.txt into all drugs
# TODO evaluate better source of features to avoid mean of means
load(here("notebooks/SCOPEAnalysis/data/well_features.RData"))
image_features <- read_csv(here("data/interim/FeatureAnalysis/feature_annotation.csv"), 
                                 col_names = TRUE)

print("preparing data")
# the overwhelming majority of samples has 1-2 probes per gene on the HG-U133 chip.
# I am aggregating intensity values using the mean to have one gene intensity per chip
expression <- promise_expr %>% 
  group_by(line, symbol, chip_name, batch, id) %>% 
  summarise(expr = mean(expr)) %>% 
  filter(!(line %in% c("D015T01", "D052T01", "D021T01")))

# I wrangle the data to fit into the MOFA input format
# TODO mean of means does not equal the overall mean across plates
morphology <- well_features %>% 
  select(one_of(paste0(image_features$name, "_expected"))) %>%
  cbind(well_metadata %>% janitor::clean_names(), .) %>% 
    as_tibble() %>% 
  dplyr::select(-plate, -well) %>%
  mutate(concentration = ifelse(drug == "DMSO", NaN, concentration)) %>%
  group_by(drug, line, concentration, replicate) %>% 
  summarise_all(mean)

# defining the set of included lines, I know there are fewer lines in the imaging data than in our gene expression corpus
lines <- morphology$line %>% unique()

# loading genetic data 
mut <-  mutations %>% 
  dplyr::select(line, symbol = SYMBOL) %>% 
  mutate(value = 1) %>% 
  complete(line, symbol, fill = list(value = 0)) %>% 
  # I discard AC008575.1, as it is a protein of unknown function overlapping with a set of exons in the APC gene. Every sample in the current dataset that has a mutation in this uncharacterized protein also has a mutation in APC
  filter(symbol != "AC008575.1") %>% 
  filter(line %in% lines) # I expect some missing data here as some lines have not been characterized thoroughly enough

df <- mut %>% filter(line == "D020T01") %>% mutate(line = "D020T02")
mut <- rbind(mut, df)

print("defining groups")
# reading and wrangling metadata from imaging to define batches
# TODO needs more rigorous joining by plate_id
morph_groups <- morph_metadata %>% 
  mutate(Line = paste0(donor, tumor)) %>% 
  filter(image_CTG == "Imaging") %>%
  dplyr::select(Line, Plate = barcode, screen_ID) %>% 
  #cleaning data to harmonize different naming schemes
  mutate(Plate = if_else(Line == "D020T01" & screen_ID %in% c("HC1092-9", "HC1092-10"), str_replace(Plate, "D020T01", "D020T02"), Plate)) %>%
  mutate(Line = if_else(Line == "D020T01" & screen_ID %in% c("HC1092-9", "HC1092-10"), "D020T02", Line)) %>% 
  dplyr::select(-Plate) %>% 
  janitor::clean_names() %>%
  distinct(line, screen_id) %>% 
  arrange(line) %>% 
  filter(line %in% lines) %>%
  cbind(replicate = rep(c(1,2), times = length(lines))) %>%
  mutate(img_group = str_split(screen_id, pattern = "-") %>% unlist() %>% tail(1) %>% as.numeric()) # used to be factor(screen_id) %>% as.numeric()
  
# preparing data
morph <- morphology %>% 
  mutate(observation = paste0(line, "_", drug, "_", concentration)) %>% 
  rename(group = replicate) %>% 
  ungroup() %>%
  dplyr::select(-drug, -concentration)
  
expr <- expression %>% 
  mutate(line = paste0(line, "T01")) %>%
  filter(line %in% lines)

expr_groups <- expr %>% 
  ungroup() %>%
  dplyr::select(line, batch) %>% 
  distinct() %>% 
  mutate(expr_group = factor(batch) %>% as.numeric()) %>% 
  rbind(tibble(line = c("D020T02", "D020T02"),
               batch = c("batch4", "batch5"),
               expr_group = c(4, 5))) 

all_groups <- morph_groups %>% 
  arrange(line) %>%
  cbind(expr_groups %>% arrange(line) %>% rename(line_expr = line)) %>% 
  mutate(overall_group = paste0(img_group, "_", expr_group)) %>% 
  dplyr::select(line, replicate, batch, overall_group, img_group, expr_group)
  
all_groups %>% 
  mutate(id = paste0(line, "_", overall_group)) %>%
  dplyr::select(id, img_group, expr_group) %>%
  as.data.frame() %>% column_to_rownames("id") %>% 
  pheatmap::pheatmap()

all_groups %>% write_rds(here("data/processed/PhenotypeSpectrum/groups.Rds"))

print("preparing final df")
# Creating a df with the following variables sample   group          feature   view value
mofa_df <- rbind(
mut %>% mutate(sample = paste0(line, "_DMSO_NaN"), feature = symbol, view = "mutation") %>% dplyr::select(sample, feature, view, value) %>% 
  mutate(group = NA),
morph %>% pivot_longer(contains("x."), names_to = "feature", values_to = "value") %>% dplyr::rename(replicate = group) %>% 
  left_join(all_groups %>% dplyr::select(line, replicate, overall_group)) %>% 
  dplyr::select(-line, -replicate) %>% 
  dplyr::rename(group = overall_group,
                sample = observation) %>% 
  mutate(view = "morphology"),
expr %>% ungroup() %>% mutate(sample = paste0(line, "_DMSO_NaN"), feature = symbol, view = "expression") %>% 
  left_join(all_groups %>% dplyr::select(line, batch, overall_group)) %>% 
  dplyr::rename(group = overall_group,
                value = expr) %>% 
  mutate(sample = paste0(line, "_DMSO_NaN")) %>% 
  dplyr::select(sample, feature, value, view, group)
)

# writing mofa_df
mofa_df %>% write_rds(here("data/processed/PhenotypeSpectrum/mofa_df.Rds"))
