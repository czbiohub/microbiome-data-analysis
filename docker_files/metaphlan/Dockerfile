# Base Image
FROM continuumio/miniconda:4.5.12

USER root:root

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Add channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install Sourmash
RUN conda install -y metaphlan2 python=2
RUN conda install -c conda-forge -y awscli python=2
RUN conda install -c bioconda -y bbmap python=2

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda:4.5.12" \
      container.version="0.0.1" \
      software.name="Metaphlan" \
      software.version="2.7.7" \
      software.description="MetaPhlAn is a computational tool for profiling the composition of microbial communities from metagenomic shotgun sequencing data with species level resolution." \
      software.website="https://bitbucket.org/biobakery/metaphlan2" \
      software.documentation="https://bitbucket.org/biobakery/metaphlan2" \
      container.category="tool"

COPY . .
RUN chmod -R +x ./