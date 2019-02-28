# Base Image
FROM continuumio/miniconda:4.5.12

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Add channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install MIDAS
RUN conda install -c bioconda -y midas vsearch hmmer
RUN conda install openssl=1.0

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda:4.5.12" \
      software.name="MIDAS" \
      software.description="An integrated pipeline for estimating strain-level genomic variation from metagenomic data" \
      software.website="https://github.com/snayfach/MIDAS" \
      software.documentation="https://github.com/snayfach/MIDAS/blob/master/README.md" \
      container.category="aligner"
      