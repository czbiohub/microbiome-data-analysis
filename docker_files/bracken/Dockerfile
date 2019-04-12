# Base Image # Python2
FROM continuumio/miniconda:4.5.12

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda:4.5.12" \
      container.version="0.0.1" \
      software.name="bracken" \
      software.version="" \
      software.description="Bracken (Bayesian Reestimation of Abundance with KrakEN) is a highly accurate statistical method that computes the abundance of species in DNA sequences from a metagenomics sample." \
      software.website="http://ccb.jhu.edu/software/bracken/index.shtml" \
      software.documentation="http://ccb.jhu.edu/software/bracken/index.shtml?t=manual" \
      container.category="species abundance"

USER root:root

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Add channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install Sourmash
RUN conda install -c conda-forge -y awscli python=2.7
RUN conda install -c bioconda -y kraken2 bracken bbmap python=2.7

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
RUN chmod -R +x ./
