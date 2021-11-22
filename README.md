# PROMISE: Image-based Profiling of Cancer Organoids

![image](https://user-images.githubusercontent.com/18559148/136934819-e9881ecf-9f37-4d30-85aa-07b2418317e6.png)

This repository contains all notebooks and supporting code for the Boutros lab cancer organoid image-based profiling project. The latest version of the corresponding manuscript is published on [bioarXiv](https://www.biorxiv.org/content/10.1101/660993v1.full). The manuscript is currently undergoing peer-review. If you have comments, criticism or questions, please feel free to create an issue in this repository, send us a [tweet](https://twitter.com/Niklas_TR) or comment at bioarXiv. 

This repository comes with a matching [docker](https://www.docker.com/products/docker-desktop) container [image](https://hub.docker.com/r/niklastr/promise/tags), which contains all dependencies and additional raw data to re-run the analysis.

The repository structure is based on the [cookiecutter datascience](https://github.com/drivendata/cookiecutter-data-science) standard and the [rocker](https://www.rocker-project.org/) docker template for R. In order to run the analysis, pull this repository from github and install the [SCOPEAnalysis](https://figshare.com/s/e465d65a9964d3b999e9) package. Alternatively, pull the pre-built docker [image](https://hub.docker.com/layers/158839806/niklastr/promise/latest/images/sha256-362bac7f1dc8bafa2bfb519413ed08ed1ec4023171cf618c17e47eca0686fbf7?context=repo) (recommended) which has the repository, the package and most dependencies preinstalled. You can interact with Rstudio Server which is running in the docker container using your [local browser](localhost:8080). Running Rstudio server from the docker image might require restarting the current R-session or starting in safe mode for some users. 

Sequencing and gene expression data has been deposited in public repositories, such as [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117548) and [EGA](https://ega-archive.org/studies/EGAS00001003140).

## Docker Containers
The project can most easily be reproduced by pulling these docker containers: 

* [niklastr/promise:latest](https://hub.docker.com/r/niklastr/promise/tags) - contains the promise git project together with large local files stored under *localdata*
* [niklastr/MOFA:latest](https://hub.docker.com/r/niklastr/mofa/tags) - contains a MOFA2 implementation to run the multi-omics factor analysis. The code can be run without GPU support
* [niklastr/promise:interimdata](https://hub.docker.com/r/niklastr/promise/tags) - contains all contents of the niklastr/promise:latest with additional interim data from the original analysis (available Q4 2021)

## Reproducing Notebooks
### Figures
All notebooks can be reproduced by starting a niklastr/promise docker, navigating to the promise working directory and calling

```
Rscript --vanilla /home/rstudio/promise/make_results.R
```

Knitted vignettes will appear in the notebook subdirectories. Individual figures are exported into the reports/figures directory.

### MOFA modeling
The MOFA model can be trained using a niklastr/MOFA docker. 
In order to run the model, the promise git directory will have to be mounted as a volume, for example via **docker run -d -v /Users/myuser/github/promise:/promise niklastr/mofa**. Once running, call the below command to initiate the data preparation and MOFA modeling. 

```
Rscript --vanilla /promise/src/models/mofa/tidy_mofa.R
```
The model will be trained and a **model.hdf5** will be created under **promise/models/mofa**.

## Directory Structure

```
├── LICENSE
├── Makefile           <- Makefile with commands like `make data` or `make train`
├── README.md          <- The top-level README for developers using this project.
├── data
│   ├── external       <- Data from third party sources including in-house R packages.
│   ├── interim        <- Intermediate data that has been transformed.
│   ├── processed      <- The final, canonical data sets for modeling.
│   └── raw            <- The original, immutable data dump. Data in this directory is incomplete due to size constraints.
│
├── models             <- Trained models, including MOFA models and UMAP clustering results
│
├── notebooks          <- R markdown notebooks. Naming convention is a number (for ordering),
│                         the creator's initials, and a short `-` delimited description, e.g.
│                         `1.0-jqp-initial-data-exploration`. Notebooks are sorted by topic.
│
├── references         <- Data dictionaries, manuals, and all other explanatory materials such as supplementary tables.
│
├── reports            <- Generated graphics and figures to be used in reporting
│
├── results            <- Camera-ready figures as pdf and illustrator files
│
├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
│                         generated with `pip freeze > requirements.txt`
│
├── src                <- Source code for use in this project.
│   │
│   ├── data           <- Scripts to download, generate and pre-process data
│   │
│   ├── models         <- Scripts to extract features, train models, create UMAP embeddings and run MOFA        
│   │
│   └── visualization  <- Scripts to create exploratory visualizations (only works with raw data access)
│
└── docker            <- dockerfiles for niklastr/promise:clean and niklastr/mofa:clean
```

## Directory Details 
### data
#### data/raw/PROMISE (non-public imaging data)
* hdf5dir - contains raw well-level hdf5 imaging data of each well (e.g. per object: 4 regions, 3 channels, 16 z-levels), ca. 1GB/well
* hdf5projection - contains hdf5 objects with image projection information on a well level, ca. 50MB/well
* hdf5validation - text files containing error reports in case checksum are not matching between hdf5 files after file transfer (data was originally moved after local projection to a remote compute environment and storage)
* htmldir - visualization of plate projections
* segmentation - image segmentation segmentation masks on well-level
* segmentationhtmldir - visualization of segmentation masks
* features - storage of feature extraction results on a plate level with 3 elements
	* wells/ - directory with well-level feature data
	* features - a plate level aggregation of raw well-level features
	* processedFeatures - these data have been pre-processed, including the following steps (see below, well-level aggregates exist in _averaged_features)
		* removing objects at the well boundary
		* removing objects that are out of focus
		* removing objects with unexpected size

#### data/interim (public data, **line_differences**  and **drug_effects** will be available via docker niklastr/promise:interimdata)
the directory contains 3 larger projects that were started for particular purposes within the manuscript: 
* **line_differences** - a collection of common features across all organoid lines, followed by PCA for Unsupervised Learning and EDA
	* results/ReducedFeatures_all_drugs_human.h5 - 27GB large file containing all shared features across lines
	* results/ReducedFeaturesPCA_all_drugs_human.h5 - 5G large file containing results of a PCA with 25 preserved principal components across all analyzed lines
* **drug_effects** - a project to estimate the effect of each drug treatment on every individual organoid line. Common features across lines are identified, features for each line are scaled to their DMSO control and PCA is being performed.
	* TransformedFeatures_human_25components.h5 - 5G large file containing results of a incremental PCA with 25 preserved components across all analyzed lines, the **important difference to ReducedFeaturesPCA_all_drugs_human.h5** is that input features were scaled to each line's respective DMSO control.
	* incrementalPCA - pickle of PCA model
	* lines/ - containing the logistic regression models for each line and drug with their respective AUROC and pvalue
* **organoid_viability**
	* classifiers/ - model checkpoint of random forest classifiers
	* diagnostics/ 
	* results/ - containing csv files with viability estimates

#### data/processed (public data, available via docker niklastr/promise:latest and git repository)
* umap_ ..
	* umap_absolute_all_drugs_tidy.Rds - ca. 1GB large file with all included organoid objects as UMAP projection with metadata
	* umap_absolute_all_drugs_sampled.Rds - a subset of the UMAP object above, that represents 5% randomly sampled organoids of the original corpus, ca. 300k objects


### src - source code required to process data and run models
#### src/data
the directory contains two libraries (PROMISE, S3O and SCOPEAnalysis). The first two packages are required for the generation of (image and feature-level) raw data. 
* PROMISE - contains a scheduler and handles the processing of raw microscopy images, incl. the compression into hdf5 objects
* S3O - is a segmentation toolbox for organoid projections and, together with PROMISE is required to generate the feature-level raw data
* make_expression.R - R script to perform preprocessing on microarray data, starting with publicly deposited .CEL files
#### src/models
code within the models directory is used to process feature-level raw data (input in data/raw, output in data/interim and data/processed)
* FeatureAnalysis - the library contains python scripts to process the raw feature data by filtering and scaling it
  * Utils.py - contains two functions that contain the manual annotation for every imaged plate and need adjustment when running on new data 
    * get_rep_ids
    * get_rep_ids_real
  * PCAScaledFeatures.py - a set of function required for the PCA projection of raw feature-level organoid data; the results are stored within the **line_differences** subdir (data/interim/line_differences)
* bsub - contains .bsub jobs for the DKFZ ODCF cluster, the scripts call other functions that can be run in a different compute environemnt as well
* umap - contains code for the further embedding of PCA projections from (ReducedFeaturesPCA_all_drugs_human.h5); the final result are the umap files **umap_absolute_all_drugs_tidy** and the downsampled version **umap_absolute_all_drugs_sampled**
* clustering contains code related to various clustering methods for UMAP or PCA level feature data
* mofa - contains code for data preparation and modeling. The central script is **tidy_mofa.R**
