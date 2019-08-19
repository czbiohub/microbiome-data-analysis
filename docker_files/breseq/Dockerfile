# Base Image
FROM continuumio/miniconda3:4.6.14 

USER root:root

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Install Bowtie2
RUN conda install -c bioconda -y breseq
RUN conda install -c conda-forge -y awscli
RUN conda clean -ya

# Script location
RUN mkdir -p /mnt /data
WORKDIR /mnt

COPY . .
RUN chmod -R +rx ./

# Data location
WORKDIR /data

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.6.14" \
      container.version="0.0.1" \
      software.name="BreSeq" \
      software.description="A computational pipeline for finding mutations relative to a reference sequence in short-read DNA re-sequencing data." \
      software.website="https://github.com/barricklab/breseq/tree/v0.33.2" \
      software.documentation="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml" \
      container.category="tool"
      