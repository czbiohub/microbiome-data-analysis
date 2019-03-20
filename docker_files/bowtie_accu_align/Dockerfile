# Setup base image as dev
FROM sunitjain/bowtie2:latest as dev

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

# Install dependencies
RUN conda install scipy numpy pandas awscli
RUN conda install -c bioconda -y samtools bedtools vcftools pysam sambamba picard fastqc bwa bbmap
RUN conda install openssl=1.0

RUN mkdir -p /mnt
WORKDIR /mnt

# Get Repo
COPY . .

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/bowtie2" \
      software.name="accu_align" \
      software.description="Brian Yu's alignment pipeline" \
      software.website="" \
      container.category="aligner"

RUN chmod -R +rx ./
