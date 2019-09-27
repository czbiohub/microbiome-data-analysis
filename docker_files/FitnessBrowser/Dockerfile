# Setup base image
FROM continuumio/miniconda3:4.6.14 

USER root:root
ENV PATH "/opt/conda/bin:$PATH"

RUN conda install -c bioconda -y perl-bioperl
RUN conda install -c conda-forge -y awscli procps-ng time
RUN conda install -c r -y r-base 
RUN conda clean -ya && \
    mkdir -p /mnt /data

RUN R -e "install.packages(c('parallel'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

WORKDIR /mnt

# Get BLAT
RUN curl -O http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod +x blat && \
    mv blat /usr/bin/

# Get Repo
COPY . .

# Get FitnessBrowser from git https://github.com/pflu-sbw25/FitnessBrowser
# Branch: Master 
# commit: 9c85d32a4f0b9baca7e27604d20cbbc0a45bce33 
# dated: March 21, 2019
RUN tar xzf FitnessBrowser.tar.gz && \
    rm FitnessBrowser.tar.gz && \
    mv scripts/* FitnessBrowser/bin/

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="continuumio/miniconda3:4.6.14" \
      software.name="FiBo" \
      software.description="Fitness Browser: Scripts for Estimating Mutant Fitness By Sequencing Randomly Barcoded Transposons" \
      software.website="https://github.com/pflu-sbw25/FitnessBrowser" \
      container.category="aligner"

RUN chmod -R +rx ./

WORKDIR /data