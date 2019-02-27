# Base Image
FROM sunitjain/bowtie2:20190227

# File Author / Maintainer
MAINTAINER Sunit Jain <sunit.jain@czbiohub.org>

# Metadata
LABEL container.base.image="sunitjain/bowtie2:20190227"
LABEL software.name="IGGSearch"
LABEL software.description="Metagenomic species profiling with enhanced coverage of the human gut microbiome"
LABEL software.website="https://github.com/snayfach/IGGsearch"
LABEL software.documentation="https://github.com/snayfach/IGGsearch/blob/master/README.md"
LABEL tags="omics/alignments"
LABEL tags="type/analysis"

# Install dependencies
RUN conda install samtools
RUN conda install numpy
RUN conda install pysam

# Install IGGSearch
RUN git clone https://github.com/snayfach/IGGsearch
RUN chmod -R +x IGGsearch

# Set environment variables
ENV PYTHONPATH=$PYTHONPATH:${PWD}/IGGsearch/iggsearch
ENV PATH=$PATH:${PWD}/IGGsearch

# Download IGGdb
RUN mkdir databases && \
  cd databases && \
  run_iggsearch.py download --gut-only && \
  tar -zxvf iggdb_v1.0.0_gut.tar.gz && \
  cd -

ENV IGG_DB=${PWD}/databases/iggdb_v1.0.0_gut

CMD ["python", "run_iggsearch.py", "search", "--m1", "test.fastq.gz", "--outdir", "test", "--test"]
