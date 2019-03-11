# Base Image. Miniconda is build on ubuntu 16.04LTS
FROM continuumio/miniconda:4.5.12

ENV PATH "/opt/conda/bin:${PATH}"

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Add channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install awscli
RUN conda install awscli
# Install MIDAS
RUN conda install -c bioconda -y bbmap fastqc

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda:4.5.12" \
      software.name="Read QC" \
      software.description="A docker image with common Read QC tools for *omics data" \
      container.category="qc"

RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
RUN echo ${PATH}