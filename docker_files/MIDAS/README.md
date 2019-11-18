# MIDAS docker

## Minimum Requirements

```{bash}
vcpus=4
memory=8000 (Mb)
size of default DB: 36GB
```

## Expected Variables

### Required

```{bash}
coreNum=4
fastq1="s3://PATH/TO/FORWARD.fastq"
S3OUTPUTPATH="s3://PATH/TO/OUTPUT/DIR/"
```

### Optional (can be left blank)

```{bash}
fastq2="s3://PATH/TO/REVERSE.fastq" OR ""
hardTrim=80 OR ""
subsetReads=4000000 OR ""
mapid=94 OR "" (applicable to Midas SNP only)
aln_cov=0.75 OR "" (applicable to Midas SNP only)
s3path2db="s3://PATH/TO/MIDAS_DB/TARBALL" OR empty -> DEFAULT: "s3://czbiohub-brianyu/Synthetic_Community/Genome_References/MIDAS/1.2/midas_db_v1.2.tar.gz"
```

### Additional options

#### Skipping Fastq QC

To skip BBDuk QC and directly use input fastqs as inputs for Midas add the following statement to the commands flag. ***Note** that this option is not yet supported by the submission script*

```{bash}
export skip_qc=true;
```

## Sample Batch Submission

```{bash}
aegea batch submit --queue aegea_batch --image sunitjain/midas:latest \
--storage /mnt=500 --memory 8000 --vcpus 4 \
--command="export coreNum=4; \
export fastq1=s3://PATH/TO/FORWARD.fastq; \
export fastq2=s3://PATH/TO/REVERSE.fastq; \
export S3OUTPUTPATH=s3://PATH/TO/OUTPUT/DIR/; \
export subsetReads=''; \
export hardTrim=''; \
export minid=''; \
export aln_cov=''; \
./run_midas.sh"
```

## Bulk Submissions

Use the [create_submission_commands.py](https://github.com/czbiohub/microbiome-data-analysis/blob/master/docker_files/submissions/create_submission_commands.py) script along the [pipelines.json](https://github.com/czbiohub/microbiome-data-analysis/blob/master/docker_files/submissions/pipelines.json) file to generate a list of aegea submission commands based on either a `seed file` or an input `S3path`.

To execute the generated list of commands, try the following:

1. Ensure you have GNU Parallel installed (try using `conda` or `homebrew`).
2. Use the following command, assuming your generated list of aegea commands is in `submissions.sh`

```{bash}
cat submissions.sh | parallel --joblog submissions.joblog --retries 3 --delay 30s "{}" &> submissions.log &
```

This will execute each line in the `submissions.sh` file with 30 seconds delay between them. This delay is helpful to avoid AWS submission errors and selection of a reasonable instance for your jobs. You may increase or decrease this number as you please or based on the number of commands you'll be executing.
