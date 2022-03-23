FROM rocker/verse:4.0.0

RUN apt-get update && apt-get install -y \
  build-essential \
  libglpk40 \
  tmux \
  libxml2-dev \
  libcurl4-openssl-dev \
  libssl-dev \ 
  nano \
  less

# R packages
## CRAN
#RUN R -e "install.packages('tidyverse',dependencies=TRUE)" 
RUN R -e "install.packages('igraph',dependencies=TRUE)" 
### grouped packages where I do not expect errors
RUN R -e "install.packages(c('janitor', 'ggdag', 'here', 'cowplot', 'gridExtra', 'pheatmap', 'princurve', 'scico'), dependencies = TRUE)"
RUN R -e "install.packages(c('ggrastr'), dependencies = TRUE)"

## Bioconductor
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('org.Hs.eg.db',dependencies=TRUE)"
RUN R -e "BiocManager::install('biomaRt',dependencies=TRUE)"
RUN R -e "BiocManager::install('enrichplot',dependencies=TRUE)"
RUN R -e "BiocManager::install('clusterProfiler',dependencies=TRUE)"
RUN R -e "BiocManager::install('ReactomePA',dependencies=TRUE)"
### grouped packages where I do not expect errors
RUN R -e "BiocManager::install(c('lumi', 'limma', 'affy'),dependencies=TRUE)"

## devtools and dependencies
### installing dependencies for ggrastr
RUN apt-get install -y \
  libgtk2.0-dev \
  libcairo2-dev \
   xvfb \
   xauth \
   xfonts-base \
   libxt-dev
# RUN R -e "devtools::install_github('VPetukhov/ggrastr', build_vignettes = FALSE)"
# RUN R -e "install.packages('ggrastr',dependencies=TRUE)" 
### installing dependencies for MOFA2
RUN R -e "devtools::install_version('matrixStats', version = '0.60.0', repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('MOFA2',dependencies=TRUE)"

## promise specific packages
RUN R -e "BiocManager::install('Biobase',dependencies=TRUE)"
RUN R -e "BiocManager::install('GEOquery',dependencies=TRUE)"
RUN R -e "BiocManager::install('arrayQualityMetrics',dependencies=TRUE)"
RUN R -e "install.packages(c('platetools'), dependencies = TRUE)"
RUN R -e "BiocManager::install('preprocessCore', configure.args='--disable-threading')"
RUN R -e "BiocManager::install('DESeq2',dependencies=TRUE)"

RUN apt-get install -y \
  libfftw3-3
RUN R -e "BiocManager::install('EBImage',dependencies=TRUE)"
RUN wget https://bioconductor.statistik.tu-dortmund.de/packages/3.11/bioc/src/contrib/MaxContrastProjection_1.11.0.tar.gz
RUN R CMD INSTALL MaxContrastProjection_1.11.0.tar.gz
#RUN mv MaxContrastProjection_1.11.0.tar.gz data/external 
# R CMD INSTALL data/external/SCOPEMouse
COPY . /home/rstudio/promise
RUN R CMD INSTALL /home/rstudio/promise/code/data/PROMISE
WORKDIR /home/rstudio/promise
# TODO fix commented PROMISE package install
# RUN Rscript --vanilla /home/rstudio/promise/make_results.R
# install ggrastr in the most complicated way
RUN wget https://cran.r-project.org/src/contrib/ragg_1.2.2.tar.gz
RUN wget https://cran.r-project.org/src/contrib/ggrastr_1.0.1.tar.gz
RUN wget https://cran.r-project.org/src/contrib/vipor_0.4.5.tar.gz
RUN wget https://cran.r-project.org/src/contrib/beeswarm_0.4.0.tar.gz
RUN wget https://cran.r-project.org/src/contrib/ggbeeswarm_0.6.0.tar.gz
RUN wget https://cran.r-project.org/src/contrib/textshaping_0.3.6.tar.gz
RUN wget https://cran.r-project.org/src/contrib/cpp11_0.4.2.tar.gz
RUN wget https://cran.r-project.org/src/contrib/systemfonts_1.0.4.tar.gz

RUN R CMD INSTALL systemfonts_1.0.4.tar.gz
RUN R CMD INSTALL vipor_0.4.5.tar.gz
RUN R CMD INSTALL beeswarm_0.4.0.tar.gz
RUN R CMD INSTALL ggbeeswarm_0.6.0.tar.gz
RUN R CMD INSTALL cpp11_0.4.2.tar.gz 
RUN R CMD INSTALL textshaping_0.3.6.tar.gz
RUN R CMD INSTALL ragg_1.2.2.tar.gz
RUN R CMD INSTALL ggrastr_1.0.1.tar.gz
CMD [ "bash" ]


