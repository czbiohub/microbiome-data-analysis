#!/bin/bash -x

set -e
set -u
set -o pipefail

TASK="cd-hit-est"
THREADS=15
DOCKER_IMAGE=quay.io/biocontainers/cd-hit
DOCKER_VERSION=4.8.1--h8b12597_3

QUERY_PATH=${1}
CLUST_ID=${2}
SHORT_SEQ_ALN=${3}

RESULTS_PATH=$(echo "${QUERY_PATH}" | sed -e 's/queries/results/' -e 's/.fasta//')

BASE_PATH=$(pwd)

LOCAL_QUERY_FASTA="${BASE_PATH}/${QUERY_PATH}"
LOCAL_RESULT_FILE="${BASE_PATH}/${RESULTS_PATH}"

LOCAL_QUERY_DIR="$(dirname ${LOCAL_QUERY_FASTA})"
QUERY_FILE_NAME="$(basename ${LOCAL_QUERY_FASTA})"

LOCAL_RESULTS_DIR="$(dirname ${LOCAL_RESULT_FILE})"
RESULTS_FILE_NAME="$(basename ${LOCAL_RESULT_FILE})__c${CLUST_ID}__aS${SHORT_SEQ_ALN}.rep.fasta"

mkdir -p ${LOCAL_RESULTS_DIR}
DOCKER_QUERIES_DIR="/cdhit/queries"
DOCKER_RESULTS_DIR="/cdhit/results"

docker container run --rm \
    -v ${LOCAL_QUERY_DIR}:${DOCKER_QUERIES_DIR}:ro \
    -v ${LOCAL_RESULTS_DIR}:${DOCKER_RESULTS_DIR}:rw \
    ${DOCKER_IMAGE}:${DOCKER_VERSION} \
        ${TASK} \
            -T ${THREADS} \
            -M 0 \
            -d 0 \
            -c ${CLUST_ID} \
            -s 0.9 \
            -aS ${SHORT_SEQ_ALN} \
            -i ${DOCKER_QUERIES_DIR}/${QUERY_FILE_NAME} \
            -o ${DOCKER_RESULTS_DIR}/${RESULTS_FILE_NAME}

# -c 0.99 \ # cluster id 99%
# -s 0.9 \ # shorter seq in a cluster must be at least 90% of longest
# -aS 0.99 \ # 99% of shorter seq must be aligned to longest
