# Base Image
FROM sunitjain/bowtie2:20190227

# Install dependencies
RUN conda install numpy pysam
RUN conda install -c bioconda -y samtools
RUN conda install openssl=1.0

# Data location
RUN mkdir /app
WORKDIR /app

# Install IGGSearch
RUN git clone https://github.com/snayfach/IGGsearch && \
    chmod -R +x IGGsearch

# Set environment variables
ENV PYTHONPATH=$PYTHONPATH:/app/IGGsearch/iggsearch
ENV PATH=$PATH:/app/IGGsearch

# Download IGGdb
RUN mkdir databases && \
    cd databases && \
    run_iggsearch.py download --gut-only && \
    tar -zxvf iggdb_v1.0.0_gut.tar.gz && \
    cd -

ENV IGG_DB=/app/databases/iggdb_v1.0.0_gut

# Metadata
LABEL container.maintainer="Sunit Jain <sunit.jain@czbiohub.org>" \
      container.base.image="sunitjain/bowtie2:20190227" \
      software.name="IGGsearch" \
      software.description="Metagenomic species profiling with enhanced coverage of the human gut microbiome" \
      software.website="https://github.com/snayfach/IGGsearch" \
      software.documentation="https://github.com/snayfach/IGGsearch/blob/master/README.md" \
      container.category="aligner"

CMD python /app/IGGsearch/run_iggsearch.py search --m1 /app/IGGsearch/test.fastq.gz --outdir /app/test --test
