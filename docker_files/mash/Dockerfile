# Base Image
FROM continuumio/miniconda3

USER root:root

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Install Sourmash
RUN conda install -c conda-forge -y awscli python=3
RUN conda install -c bioconda -y mash python=3

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.5.12" \
      container.version="0.0.1" \
      software.name="Mash" \
      software.version="2.0.0a9" \
      software.description="Compute and compare MinHash signatures for DNA data sets." \
      software.website="https://github.com/marbl/Mash" \
      software.documentation="http://mash.readthedocs.org" \
      container.category="kmer sketch"

COPY . .
RUN chmod -R +x ./