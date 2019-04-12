# Base Image
FROM sunitjain/bowtie2:latest

USER root

# Install dependencies
RUN conda install numpy
RUN conda install -c bioconda -y samtools pysam bbmap
RUN conda install -c conda-forge -y awscli
RUN conda install openssl=1.0

RUN mkdir -p /mnt

# map /mnt to /efs
# VOLUME "/mnt/databases"

# Data location
WORKDIR /mnt

# Install IGGSearch
# RUN git clone --depth 1 https://github.com/snayfach/IGGsearch && \
#     chmod -R +x IGGsearch

# No longer required if mounting the EFS or downloading from S3
# Download IGGdb
# RUN cd databases && \
#     aws s3 cp s3://czbiohub-microbiome/IGGdb/gut_only/v1.0.0/iggdb_v1.0.0_gut.tar.gz . && \
#     tar xzf iggdb_v1.0.0_gut.tar.gz && \
#     cd -

COPY . .
RUN chmod -R +rx ./

# Set environment variables
# ENV PYTHONPATH "/mnt/bin/iggsearch"
# ENV PATH "/mnt/bin:/opt/conda/bin:${PATH}"
# ENV IGG_DB "/mnt/databases/iggdb_v1.0.0_gut"

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="sunitjain/bowtie2:20190227" \
      software.name="IGGsearch" \
      software.description="Metagenomic species profiling with enhanced coverage of the human gut microbiome" \
      software.website="https://github.com/snayfach/IGGsearch" \
      software.documentation="https://github.com/snayfach/IGGsearch/blob/master/README.md" \
      container.category="aligner"
