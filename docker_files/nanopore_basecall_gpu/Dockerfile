# Base Image
FROM nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04
USER root:root

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

# Install system-level programs
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install --no-install-recommends --no-install-suggests -qqy apt-transport-https wget tzdata bzip2 ca-certificates curl git && \
    apt-get install --no-install-recommends --no-install-suggests -qqy lsb-release pigz gnupg gnupg2 software-properties-common apt-utils

# Install Oxford Nanopore software
RUN export PLATFORM=$(lsb_release -cs) && \
    wget -O- https://americas.oxfordnanoportal.com/apt/ont-repo.pub | apt-key add - && \
    echo "deb http://americas.oxfordnanoportal.com/apt ${PLATFORM}-stable non-free" | tee /etc/apt/sources.list.d/nanoporetech.sources.list && \
    apt-get update && \
    apt-get install --no-install-recommends --no-install-suggests -qqy ont-guppy libcairo2-dev nvidia-384- libcuda1-384- && \
    apt-get update --fix-missing && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Installing miniconda3, copied from https://hub.docker.com/r/continuumio/miniconda3/dockerfile
# ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    echo "export PATH=/opt/conda/bin:$PATH" >> ~/.bashrc

COPY . .

# Configuring conda channels
RUN conda init bash && \
    conda config --add channels bioconda && \
    conda update -n base -c defaults -y -q conda && \
    conda install -y -c conda-forge awscli && \
    conda install -y -c r r-base && \
    Rscript install_packages.R && \
    wget https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R -O /mnt/MinIONQC.R && \
    chown -R root * && chgrp -R root * && chmod -R +rx ./
    # Rscript -e "install.packages('futile.logger', repos='https://cran.rstudio.com/', dep=TRUE)" && \

# Metadata
LABEL container.maintainer="Brian Yu <brian.yu@czbiohub.org>" \
      container.base.image="nvidia/cuda:9.0-cudnn7-devel-ubuntu16.04" \
      container.name="nanopore_basecall" \
      container.version="0.1.1" \
      software.name="Miniconda 4.5.11, ont-guppy, R, MinIONQC"
