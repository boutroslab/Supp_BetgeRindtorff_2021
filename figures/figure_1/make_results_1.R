# run this script from within the docker container by calling Rscript --vanilla make_results.R
# the git directory and accompanying docker container come "batteries included". This script can be run to recreate all figures in one go.

# defining lib path
#.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/4.0")
print(.libPaths())

# figure 1
# requires running the line_difference scripts within the FeatureAnalysis module
# the arguments supplied within the knit command override the defaults defined within each vignette. 
# This can be valuable, for example, for rerunning the UMAP plots with the non-downsampled complete dataset 
rmarkdown::render(here::here('notebooks/imaging/1.0-nr-organoid_unsupervised_exploration.Rmd'), 
    params = list(
        remote = FALSE, 
        data = "data/processed/morphology/umap_absolute_all_drugs_sampled.Rds", 
        sample = "data/processed/morphology/umap_absolute_all_drugs_tidy_Paclitaxel.Rds",
        cache = TRUE))
