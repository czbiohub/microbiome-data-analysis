# Setup base image
FROM continuumio/miniconda

USER root:root

# Update conda to latest version.
RUN apt-get update && \
    apt-get install build-essential && \
    apt autoremove && \
    apt clean

RUN conda update -n base -c defaults conda

# Install Bowtie2
RUN conda install -y numpy scipy
RUN conda install -c bioconda -y bwa samtools bbmap
RUN conda install -c conda-forge -y awscli
RUN pip install openopt==0.5625

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
RUN chmod -R +rx ./

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda" \
      software.name="StrainFinder" \
      software.description="Find maximum likelihood estimates for the strain genotypes and the strain frequencies across all metagenomes" \
      software.website="https://github.com/cssmillie/StrainFinder" \
      container.category="read_abundance"