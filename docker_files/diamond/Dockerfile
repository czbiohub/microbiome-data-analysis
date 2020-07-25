FROM continuumio/miniconda3
ENV VERSION 2.13 

RUN apt update && \
      apt install -y procps && \
      apt-get clean && \
      rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels default

RUN conda install metabat2 && conda clean -a

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
RUN chmod -R +rx ./

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.6.14" \
      container.version="0.8.24" \
      software.name="Diamond" \
      software.description="An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences" \
      software.website="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" \
      software.documentation="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml" \
      container.category="tool"
      