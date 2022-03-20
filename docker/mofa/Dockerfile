FROM gtca/mofa2:latest

RUN apt-get update && apt-get install -y \
  build-essential \
  libglpk40 \
  tmux \
  libxml2-dev \
  libcurl4-openssl-dev \
  libssl-dev

RUN R -e "install.packages('tidyverse',dependencies=TRUE)"
RUN R -e "install.packages('here',dependencies=TRUE)"

# for promise package: install SCOPEAnalysis package
# R CMD INSTALL data/external/SCOPEAnalysis_0.1.0.tar.gz 