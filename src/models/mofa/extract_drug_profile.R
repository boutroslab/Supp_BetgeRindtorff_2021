# defining lib path
.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/3.6")
print(.libPaths())

library(rhdf5) # requried
library(tidyr)
library(dplyr)
library(magrittr)
library(here)

extract_drug_profile = function(path, drug){
    print(path)
    print('extracting average profile vector')

    h5ls(path)

    df <- h5read(path, "features_organoids") %>% as.data.frame() %>% 
    magrittr::set_colnames(h5read(path, "feature_names_organoids"))

    print(head(df))
    print(nrow(df))
}

# testing
path_in = here::here("data/interim/FeatureAnalysis/line_differences/human/all_drugs/results/filtered_lenient/ReducedFeatures_all_drugs_human.h5")

extract_drug_profile(path = path_in)
