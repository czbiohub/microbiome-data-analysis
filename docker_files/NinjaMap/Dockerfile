# Setup base image
FROM continuumio/miniconda3:4.6.14 

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

# Install dependencies
RUN conda install -y scipy numpy pandas
RUN conda install -c conda-forge -y awscli
RUN conda install -c bioconda -y bowtie2 samtools bedtools vcftools sambamba pysam bbmap
RUN conda clean -ya

RUN mkdir -p /mnt
WORKDIR /mnt

# Get Repo
COPY . .

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.6.14" \
      software.name="ninjaMap" \
      software.description="Strain abundance pipeline" \
      software.website="" \
      container.category="aligner"

RUN chmod -R +rx ./
