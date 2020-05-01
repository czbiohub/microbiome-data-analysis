#!/bin/bash -x

set -e
set -u
set -o pipefail

THREADS=15
DOCKER_IMAGE=metabat/metabat
DOCKER_VERSION=latest

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
        runMetaBat.sh \
            test/contigs.fa  \
            test/contigs-1000.fastq.bam
