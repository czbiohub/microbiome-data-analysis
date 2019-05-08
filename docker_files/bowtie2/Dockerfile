# Base Image
FROM continuumio/miniconda3:4.5.12

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Install Bowtie2
RUN conda install -c bioconda -y bowtie2 samtools bbmap
RUN conda install -c conda-forge -y awscli

# Data location
RUN mkdir -p /mnt
WORKDIR /mnt

COPY . .
RUN chmod -R +rx ./

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.5.12" \
      container.version="0.0.1" \
      software.name="Bowtie2" \
      software.description="An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences" \
      software.website="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" \
      software.documentation="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml" \
      container.category="tool"
      