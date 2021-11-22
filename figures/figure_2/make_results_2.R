# run this script from within the docker container by calling Rscript --vanilla make_results.R
# the git directory and accompanying docker container comer "batteries included". This script can be run to recreate all figures in one go.

# defining lib path
#.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/4.0")
print(.libPaths())

# figure 2
# requires running the line_difference scripts within the FeatureAnalysis module
# the arguments supplied within the knit command override the defaults defined within each vignette. 
# This can be valuable, for example, for rerunning the UMAP plots with the non-downsampled complete dataset 
# requires running the line_difference scripts within the FeatureAnalysis module
rmarkdown::render(here::here('notebooks/imaging/2.0-nr-embedding_inspection.Rmd'), 
    params = list(
        data = "data/processed/morphology/umap_absolute_all_drugs_sampled.Rds",
        data_harmony = "data/processed/morphology/harmony_umap_absolute_all_drugs_sampled.Rds",
        remote = FALSE,
        cache = TRUE))
rmarkdown::render(here::here('notebooks/drug_activity/2.0-js-OrganoidViability.Rmd'))
