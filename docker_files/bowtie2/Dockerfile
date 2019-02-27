# Base Image
FROM continuumio/miniconda3:4.5.12

# File Author / Maintainer
MAINTAINER Sunit Jain <sunit.jain@czbiohub.org>

# Metadata
LABEL container.base.image="continuumio/miniconda3:4.5.12"
LABEL software.name="Bowtie2"
LABEL software.version="2.3.4.3-1"
LABEL software.description="An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences"
LABEL software.website="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml"
LABEL software.documentation="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml"
LABEL tags="omics/alignments"
LABEL tags="type/utility"

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Add channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install awscli
RUN conda install awscli

# Install Bowtie2
RUN conda install bowtie2

CMD ["bowtie2", "--version"]
