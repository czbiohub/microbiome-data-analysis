# Setup base image as dev
FROM sunitjain/bowtie2:latest as dev

# Install dependencies
RUN conda install scipy numpy pandas
RUN conda install -c bioconda -y samtools bedtools vcftools pysam sambamba picard fastqc bwa bbmap
RUN conda install openssl=1.0

RUN mkdir -p /mnt
WORKDIR /mnt

# Get Repo
COPY . .
# coreNum=4; numPerCore=1G; 

RUN chmod -R +x ./

CMD ["/mnt/accu_align.sh"]