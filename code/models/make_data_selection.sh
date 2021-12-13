# run with R/4.0.0

# requried input for MOFA 
Rscript --vanilla src/models/mofa/aggregate_data_pca.R

# generate subsamples for vizualization
Rscript --vanilla src/models/umap/select_data_umap.R Paclitaxel 0.3 123
Rscript --vanilla src/models/umap/select_data_umap.R Bortezomib 0.3 123
Rscript --vanilla src/models/umap/select_data_umap.R Binimetinib 0.3 123
Rscript --vanilla src/models/umap/select_data_umap.R Trametinib 0.3 123
Rscript --vanilla src/models/umap/select_data_umap.R WYE-132 0.3 123
Rscript --vanilla src/models/umap/select_data_umap.R CHIR-98014 0.3 123
Rscript --vanilla src/models/umap/select_data_umap.R "Irinotecan / SN-38" 0.3 123

