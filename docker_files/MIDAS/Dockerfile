# Base Image
FROM continuumio/miniconda:4.5.12

# File Author / Maintainer
MAINTAINER Sunit Jain <sunit.jain@czbiohub.org>

# Metadata
LABEL container.base.image="continuumio/miniconda:4.5.12"
LABEL software.name="MIDAS"
LABEL software.description="An integrated pipeline for estimating strain-level genomic variation from metagenomic data"
LABEL software.website="https://github.com/snayfach/MIDAS"
LABEL software.documentation="https://github.com/snayfach/MIDAS/blob/master/README.md"
LABEL tags="omics/alignments"
LABEL tags="type/analysis"

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Add channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install dependencies
RUN conda install awscli

# Install MIDAS
RUN conda install midas

CMD ["midas"]
