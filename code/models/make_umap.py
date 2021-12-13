import os
######################################
# ONLY runs on system with IBM LSF Cluster, e.g. ssh rindtorf@odcf-lsf01.dkfz.de
# cd /dkfz/groups/shared/OE0049/B110-Isilon2/promise
print('make sure you are running make_dataset on an instance that supports IBM bsub')

# 2. create UMAP embedding on all organoid data
# ====================================

# # generate umap embedding DMSO-treated organoid features in PCA space
# os.system('bsub -R "rusage[mem=100GB]" -q long ./src/models/bsub/run_umap_dmso.bsub \
# 	-o src/models/bsub/out.bsub \
# 	-e src/models/bsub/error.bsub')
	
# # hyperparameter exploration for UMAP
# os.system('bsub -R "rusage[mem=200GB]" -q verylong ./src/models/bsub/run_umap_all_drugs_paramsearch.bsub \
# 	-o src/models/bsub/out.bsub \
# 	-e src/models/bsub/error.bsub')

# generate umap embedding with or without harmony correction on all PCA transformed organoid features
os.system('bsub -R "rusage[mem=200GB]" -q verylong ./src/models/bsub/run_umap_all_drugs.bsub \
	-o src/models/bsub/out.bsub \
	-e src/models/bsub/error.bsub')

