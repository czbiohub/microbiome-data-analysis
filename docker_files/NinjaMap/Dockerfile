# Setup base image
FROM continuumio/miniconda3:4.5.12

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

# Install dependencies
RUN conda install scipy numpy pandas
RUN conda install -c conda-forge awscli
RUN conda install -c bioconda -y bowtie2 samtools bedtools vcftools pysam bbmap

RUN mkdir -p /mnt
WORKDIR /mnt

# Get Repo
COPY . .

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.5.12" \
      software.name="accu_align" \
      software.description="Brian Yu's alignment pipeline" \
      software.website="" \
      container.category="aligner"

RUN chmod -R +rx ./
