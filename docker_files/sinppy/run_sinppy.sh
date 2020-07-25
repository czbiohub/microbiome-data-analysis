#!/bin/bash -x
# shellcheck disable=SC2086

set -e
set -u
set -o pipefail

CONTAINER_NAME="quay.io/biocontainers/snippy"
CONTAINER_VERSION="4.6.0--0"

THREADS=15

S3_REF_GENOME=${1}
S3_FWD=${2}
S3_REV=${3}
S3_OUTPUT_DIR=${4%/}

#Create the output directory
S3_GENOME_BASE=$(dirname ${S3_OUTPUT_DIR})
LOCAL_REF_BASE=$(basename ${S3_GENOME_BASE})
LOCAL_SAMPLE_BASE=$(basename ${S3_OUTPUT_DIR})

LOCAL_OUTPUT_BASE=${LOCAL_REF_BASE}/${LOCAL_SAMPLE_BASE}
RGID="${LOCAL_REF_BASE}__${LOCAL_SAMPLE_BASE}"

LOCAL_TMP="${LOCAL_OUTPUT_BASE}/TMP_${RGID}"
LOCAL_OUTPUT="${LOCAL_OUTPUT_BASE}/Snippy"
mkdir -p "${LOCAL_TMP}"

# Download
aws s3 cp ${S3_REF_GENOME} ${LOCAL_TMP}/
aws s3 cp ${S3_FWD} ${LOCAL_TMP}/
aws s3 cp ${S3_REV} ${LOCAL_TMP}/

LOCAL_REF_GENOME=$(basename ${S3_REF_GENOME})
LOCAL_FWD=$(basename ${S3_FWD})
LOCAL_REV=$(basename ${S3_REV})

# Run
docker container run --rm \
    --workdir $(pwd)/${LOCAL_OUTPUT_BASE} \
    --volume $(pwd)/${LOCAL_OUTPUT_BASE}:/output \
    --volume $(pwd)/${LOCAL_TMP}:/input \
    ${CONTAINER_NAME}:${CONTAINER_VERSION} \
    snippy \
        --cpus ${THREADS} \
        --report \
        --rgid ${RGID} \
        --outdir /output/Snippy \
        --basequal 15 \
        --ref /input/${LOCAL_REF_GENOME} \
        --R1 /input/${LOCAL_FWD} \
        --R2 /input/${LOCAL_REV}

# Cleanup downloaded files
rm -rf ${LOCAL_TMP}

# Sync Ouptut to S3
aws s3 sync ${LOCAL_OUTPUT} ${S3_OUTPUT_DIR}/Snippy