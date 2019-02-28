# Setup base image as dev
FROM sunitjain/bowtie2:latest as dev

# Install dependencies
RUN conda install scipy numpy pandas
RUN conda install -c bioconda -y samtools bedtools vcftools pysam sambamba picard fastqc bwa bbmap
RUN conda install openssl=1.0

RUN mkdir /work
WORKDIR /work

# Get Repo
RUN git clone https://github.com/czbiohub/microbiome-data-analysis.git

# coreNum=4; numPerCore=1G; 

RUN chmod -R +x ./

# ENTRYPOINT [ "bash" ]

# CMD [ "accu_align.sh", "maxInsert","maxAlignments", \
#     "genomeReferenceDir","sampleName", "trimQuality","minLength", \
#     "kmer_value","min_kmer_value", "fastq1","fastq2" ]