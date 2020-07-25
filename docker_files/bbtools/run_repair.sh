#!/bin/bash -x

set -euo pipefail

NAME=${1}
S3FWD=${2}
S3REV=${3}
S3OUTPUTPATH=${4:-"s3://czbiohub-microbiome/Sonnenburg_Lab/InfantMicrobiome/staging"}
S3OUTPUTPATH=${S3OUTPUTPATH%/}

# Setup
BBMAP_DOCKER_IMAGE="quay.io/biocontainers/bbmap"
BBMAP_IMAGE_VERSION="38.86--h1296035_0"

# Output files
FWD_OUT=${NAME}/${NAME}_R1.fastq.gz
REV_OUT=${NAME}/${NAME}_R2.fastq.gz
SIN_OUT=${NAME}/${NAME}_R0.fastq.gz

# Local fastq files
LOCAL_FWD=${NAME}/tmp/$(basename ${S3FWD})
LOCAL_REV=${NAME}/tmp/$(basename ${S3REV})

# Download files from S3
aws s3 cp ${S3FWD} ${LOCAL_FWD} --quiet
aws s3 cp ${S3REV} ${LOCAL_REV} --quiet

# Run BBTools repair.sh
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd):$(pwd)" \
    ${BBMAP_DOCKER_IMAGE}:${BBMAP_IMAGE_VERSION} \
    repair.sh \
        -Xmx20g \
        -eoom \
        in=${LOCAL_FWD} \
        in2=${LOCAL_REV} \
        out=${FWD_OUT} \
        out2=${REV_OUT} \
        outs=${SIN_OUT} \
        ziplevel=2

# Remove original files
rm -rf ${NAME}/tmp/

# Sync repaired files to S3
aws s3 sync ${NAME} ${S3OUTPUTPATH}/${NAME} --quiet

# CleanUp local directory
rm -rf ${NAME}