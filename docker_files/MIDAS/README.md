# MIDAS docker

# Requirements

  vcpus: 4
  memory: 20GB
  volume:
       default DB: 36GB

# Expected Variables

  Required:
      coreNum=4
      fastq1="s3://PATH/TO/FORWARD.fastq"
      S3OUTPUTPATH="s3://PATH/TO/OUTPUT/DIR/"

# Optional (can be left blank)

      fastq2="s3://PATH/TO/REVERSE.fastq" OR ""
      hardTrim=80 OR ""
      subsetReads=4000000 OR ""
      mapid=94 OR "" (applicable to Midas SNP only)
      aln_cov=0.75 OR "" (applicable to Midas SNP only)
      s3path2db="s3://PATH/TO/MIDAS_DB/TARBALL" OR empty -> DEFAULT: "s3://czbiohub-brianyu/Synthetic_Community/Genome_References/MIDAS/1.2/midas_db_v1.2.tar.gz"

# Sample Batch Submission

  aegea batch submit --queue aegea_batch --image sunitjain/midas:latest \
    --storage /mnt=500 --memory 8000 --vcpus 4 \
    --command="export coreNum=16; \
    export fastq1=s3://PATH/TO/FORWARD.fastq; \
    export fastq2=s3://PATH/TO/REVERSE.fastq; \
    export S3OUTPUTPATH=s3://PATH/TO/OUTPUT/DIR/; \
    export subsetReads=''; \
    export hardTrim=''; \
    export minid=''; \
    export aln_cov=''; \
    ./run_midas.sh"
