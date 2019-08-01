# Setup base image
FROM continuumio/miniconda3:4.6.14 

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

RUN conda install -c bioconda -y perl-bioperl
RUN conda install -c conda-forge -y awscli
RUN conda clean -ya && \
    mkdir -p /mnt /data

WORKDIR /mnt

# Get BLAT
RUN curl -O http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod +x blat

# Get Repo
COPY . .

# Get FitnessBrowser from git
# Branch: Master 
# commit: 9c85d32a4f0b9baca7e27604d20cbbc0a45bce33 
# dated: March 21, 2019
RUN tar xzf FitnessBrowser.tar.gz && \
    rm FitnessBrowser.tar.gz

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.6.14" \
      software.name="FiBo" \
      software.description="Fitness Browser: Scripts for Estimating Mutant Fitness By Sequencing Randomly Barcoded Transposons" \
      software.website="https://github.com/pflu-sbw25/FitnessBrowser" \
      container.category="aligner"

RUN chmod -R +rx ./

WORKDIR /data