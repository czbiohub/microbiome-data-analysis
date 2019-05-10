# Base Image
FROM ubuntu:16.04

USER root:root
# Update conda to latest version.
RUN apt-get install -y build-essential

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

COPY cellranger-3.0.2.tar.gz .

RUN tar xzf cellranger-3.0.2.tar.gz

ENV PATH=/mnt/cellranger-3.0.2:$PATH

COPY run_cellranger.sh .
RUN chmod -R +rx ./

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="ubuntu:16.04" \
      container.version="0.0.1" \
      software.name="cellRanger v3.0.2" \
      software.description="Cell Ranger is a set of analysis pipelines that process Chromium single-cell RNA-seq output to align reads, generate feature-barcode matrices and perform clustering and gene expression analysis. " \
      software.website="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger" \
      software.documentation="https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger" \
      container.category="pipeline"
      