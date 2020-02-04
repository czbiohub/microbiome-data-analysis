# Base Image
FROM nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04
USER root:root

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

# Install system-level programs
RUN apt-get update && apt-get install --no-install-recommends --no-install-suggests -qqy apt-transport-https && apt-get install --no-install-recommends --no-install-suggests -qqy wget lsb-release pigz gnupg gnupg2 software-properties-common

# Install Oxford Nanopore software
RUN export PLATFORM=$(lsb_release -cs) && wget -O- https://mirror.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add - && echo "deb http://mirror.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list && apt-get update && apt-get install --no-install-recommends --no-install-suggests -qqy ont-guppy nvidia-384- libcuda1-384-

# Installing R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/' && apt-get update && apt-get install --no-install-recommends --no-install-suggests -qqy r-base

COPY install_packages.R .

RUN Rscript install_packages.R
RUN wget https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R -O /mnt/MinIONQC.R

# Installing miniconda3, copied from https://hub.docker.com/r/continuumio/miniconda3/dockerfile
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc


# Configuring conda channels
RUN conda config --add channels defaults && conda config --add channels bioconda && conda config --add channels conda-forge

# Installing conda packages
RUN conda install -y -c conda-forge awscli

COPY . .
RUN chown -R root * && chgrp -R root * && chmod -R +rx ./

# Metadata
LABEL container.maintainer="Bryan Merrill <bmerrill@stanford.edu>" \
      container.base.image="nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04" \
      container.version="0.0.1" \
      software.name="Miniconda 4.5.11, ont-guppy, R, MinIONQC"
