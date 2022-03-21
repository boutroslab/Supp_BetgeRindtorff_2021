# PROMISE: Image-based Profiling of Cancer Organoids

![image](https://user-images.githubusercontent.com/18559148/136934819-e9881ecf-9f37-4d30-85aa-07b2418317e6.png)

## Authors
Johannes Betge1,2,3,4,\*, Niklas Rindtorff1,\*, Jan Sauer1,5,\*, Benedikt Rauscher1,\*, Clara Dingert1, Haristi Gaitantzi2,3, Frank Herweck2,3, Kauthar Srour-Mhanna2,3,4, Thilo Miersch1, Erica Valentini1, Britta Velten6, Veronika Hauber2,3, Tobias Gutting2,3, Larissa Frank1, Sebastian Belle2,3, Timo Gaiser3,7, Inga Buchholz2,3, Ralf Jesenofsky2,3, Nicolai Härtel2,3, Tianzuo Zhan1,2,3, Bernd Fischer5, Katja Breitkopf-Heinlein2,3, Elke Burgermeister2,3, Matthias P. Ebert2,3,#, Michael Boutros1,8,#
 
1 German Cancer Research Center (DKFZ), Division Signaling and Functional Genomics, and Heidelberg University, Medical Faculty Mannheim, Department of Cell and Molecular Biology, Mannheim, Germany,
2 Heidelberg University, Department of Medicine II, University Medical Center Mannheim, Medical Faculty Mannheim, Mannheim, Germany, 
3 Mannheim Cancer Center, Mannheim, Germany, 
4 German Cancer Research Center (DKFZ), Junior Clinical Cooperation Unit Translational Gastrointestinal Oncology and Preclinical Models, Heidelberg, Germany,
5 German Cancer Research Center (DKFZ), Computational Genome Biology Group, Heidelberg, Germany, 
6 German Cancer Research Center (DKFZ), Computational Genomics and Systems Genetics, Heidelberg, Germany, 
7 Heidelberg University, Institute of Pathology, University Medical Center Mannheim, Medical Faculty Mannheim, Mannheim, Germany, 
8 German Cancer Consortium (DKTK), Heidelberg, Germany
*these authors contributed equally to this study
#addresses for correspondence: m.boutros@dkfz.de and matthias.ebert@medma.uni-heidelberg.de


## Abstract
Patient derived organoids resemble the biology of tissues and tumors, enabling ex vivo modeling of human diseases from primary patient samples. Organoids can be used as models for drug discovery and are being explored to guide clinical decision making. Patient derived organoids can have heterogeneous morphologies with unclear biological causes and relationship to treatment response. Here, we used high-throughput, image-based profiling to quantify phenotypes of over 5 million individual colorectal cancer organoids after treatment with more than 500 small molecules. Integration of data using a joint multi-omics modelling framework identified organoid size and cystic vs. solid organoid architecture as axes of morphological variation observed across organoid lines. Mechanistically, we found that organoid size was linked to IGF1 receptor signaling while a cystic organoid architecture was associated with an LGR5+ stemness program. Treatment-induced organoid morphology reflected organoid viability, drug mechanism of action, and was biologically interpretable using joint modelling. Inhibition of MEK led to cystic reorganization of organoids and increased expression of LGR5, while inhibition of mTOR induced IGF1 receptor signaling. In conclusion, we identified two shared axes of variation for colorectal cancer organoid morphology, their underlying biological mechanisms and pharmacological interventions with the ability to move organoids along them. Image-based profiling of patient derived organoids coupled with multi-omics integration facilitates drug discovery by linking drug responses with underlying biological mechanisms.
Keywords:  organoids, drug profiling, high-throughput screening, cancer signaling, functional genomics, computer vision, machine learning 

## Reproducibility
This repository comes with a matching [docker](https://www.docker.com/products/docker-desktop) container [image](https://hub.docker.com/r/niklastr/promise/tags), which contains all dependencies and additional raw data to re-run the analysis.

The repository structure is based on the [cookiecutter datascience](https://github.com/drivendata/cookiecutter-data-science) standard and the [rocker](https://www.rocker-project.org/) docker template for R. In order to run the analysis, pull this repository from github and install the [SCOPEAnalysis](https://figshare.com/s/e465d65a9964d3b999e9) package. Alternatively, pull the pre-built docker [image](https://hub.docker.com/layers/158839806/niklastr/promise/latest/images/sha256-362bac7f1dc8bafa2bfb519413ed08ed1ec4023171cf618c17e47eca0686fbf7?context=repo) (recommended) which has the repository, the package and most dependencies preinstalled. You can interact with Rstudio Server which is running in the docker container using your [local browser](localhost:8080). Running Rstudio server from the docker image might require restarting the current R-session or starting in safe mode for some users. 

Sequencing and gene expression data has been deposited in public repositories, such as [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117548) and [EGA](https://ega-archive.org/studies/EGAS00001003140).

### Docker Containers
The project can most easily be reproduced by pulling these docker containers: 

* [niklastr/promise:latest](https://hub.docker.com/r/niklastr/promise/tags) - contains the promise git project together with large local files stored under *localdata*
* [niklastr/MOFA:latest](https://hub.docker.com/r/niklastr/mofa/tags) - contains a MOFA2 implementation to run the multi-omics factor analysis. The code can be run without GPU support

You can run the main docker container for the project in an Rstudio Server (does not work on M1 Mac) using the command below:

```
docker run -d  -e PASSWORD=promise -p 8080:8787 -v /Users/nrindtor/github/supp_BetgeRindtorff_2021:/home/rstudio/promise niklastr/promise:clean
```


### Reproducing Notebooks
#### Figures
All notebooks can be reproduced by starting a niklastr/promise docker, navigating to the promise working directory and calling

```
Rscript --vanilla /home/rstudio/promise/make_results.R
```

Knitted vignettes will appear in the notebook subdirectories. Individual figures are exported into the reports/figures directory.

#### MOFA modeling
The MOFA model can be trained using a niklastr/MOFA docker. 
In order to run the model, the promise git directory will have to be mounted as a volume, for example via **docker run -d -v /Users/myuser/github/promise:/promise niklastr/mofa**. Once running, call the below command to initiate the data preparation and MOFA modeling. 

```
Rscript --vanilla /promise/src/models/mofa/tidy_mofa.R
```
The model will be trained and a **model.hdf5** will be created under **promise/models/mofa**.

