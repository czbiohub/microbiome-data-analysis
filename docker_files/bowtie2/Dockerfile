# Base Image
FROM continuumio/miniconda3:4.5.12

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Add channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install Bowtie2
RUN conda install bowtie2

# Data location
RUN mkdir /app
WORKDIR /app

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.5.12" \
      container.version="0.0.1" \
      software.name="Bowtie2" \
      software.description="An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences" \
      software.website="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" \
      software.documentation="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml" \
      container.category="tool"

CMD bowtie2 --version
