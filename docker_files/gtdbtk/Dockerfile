# Setup base image
FROM continuumio/miniconda

USER root:root
RUN mkdir -p /mnt /data
WORKDIR /mnt

RUN conda install -c bioconda -y gtdbtk
RUN conda install -c conda-forge -y awscli

COPY . .
RUN chmod -R +rx ./

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda" \
      software.name="GTDBtk" \
      software.description="Genome Taxonomy Database toolkit" \
      software.website="https://gtdb.ecogenomic.org/" \
      container.category="taxonomy"