FROM continuumio/miniconda:latest

ARG DEBIAN_FRONTEND="noninteractive"

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

# Install dependencies
RUN conda install --freeze-installed -c conda-forge -yq awscli glpk \
      && conda install --freeze-installed -c bioconda -yq humann2 \
      && conda clean -afy

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
RUN chmod -R +rx ./

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda" \
      container.version="20200616" \
      software.name="Humann2" \
      software.description="HMP Unified Metabolic Analysis Network v2" \
      software.website="https://huttenhower.sph.harvard.edu/humann" \
      software.documentation="https://github.com/biobakery/humann" \
      container.category="tool"
