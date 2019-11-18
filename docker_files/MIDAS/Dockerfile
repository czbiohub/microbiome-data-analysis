# Base Image. Miniconda is build on ubuntu 16.04LTS
FROM continuumio/miniconda:4.7.12

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
RUN conda install -c bioconda -y midas vsearch hmmer bbmap
RUN conda install openssl=1.0

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda:4.7.12" \
      software.name="MIDAS" \
      software.description="An integrated pipeline for estimating strain-level genomic variation from metagenomic data" \
      software.website="https://github.com/snayfach/MIDAS" \
      software.documentation="https://github.com/snayfach/MIDAS/blob/master/README.md" \
      container.category="aligner"

USER root:root

RUN mkdir -p /mnt
WORKDIR /mnt

RUN ln -s /opt/conda/bin/hs-blastn . && \
    ln -s /opt/conda/bin/bowtie2-build . && \
    ln -s /opt/conda/bin/bowtie2 . && \
    ln -s /opt/conda/bin/samtools .

COPY . .
RUN chmod +x run_midas.sh && \
    echo ${PATH}