# run this script from within the docker container by calling Rscript --vanilla make_results.R
# the git directory and accompanying docker container comer "batteries included". This script can be run to recreate all figures in one go.

# defining lib path
#.libPaths("/omics/groups/OE0049/b110_data/B110-Isilon2/promise/x86_64-pc-linux-gnu-library/4.0")
print(.libPaths())

# # figure 3
# requires running the drug_effect scripts within the FeatureAnalysis module
rmarkdown::render(here::here('notebooks/drug_activity/3.0-js-DrugPhenotypes.Rmd'))

# figure 4, 5 and 6
# running these vignettes requires a trained MOFA model. A MOFA model can be trained from within the MOFA2 docker container
# by calling the tidy_mofa.R script
rmarkdown::render(here::here('notebooks/mofa/1.0-nr-mofa_exploration.Rmd'))
rmarkdown::render(here::here('notebooks/mofa/2.0-nr-mofa_drug_effect.Rmd'))
