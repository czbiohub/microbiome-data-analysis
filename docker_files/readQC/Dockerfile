# Base Image
FROM continuumio/miniconda3:latest

USER root:root

# Install dependencies
RUN conda install -c conda-forge -y awscli
RUN conda install -c bioconda -y bbmap fastqc
RUN conda clean -ya

RUN mkdir -p /mnt

# Data location
WORKDIR /mnt

COPY . .
RUN chmod -R +rx ./

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:latest" \
      software.name="readQC" \
      software.description="Quality Control for Illumina Reads" \
      software.website="" \
      software.documentation="" \
      container.category="QC"
