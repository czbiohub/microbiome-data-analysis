#!/bin/bash -x

set -euo pipefail

TASK=${TASK:-"blastx"}
QUERY_EXT=${QUERY_EXT:-".fasta"}
THREADS=${THREADS:-50}

QUERY_PATH=${1}
LOCAL_DB=${2:-"/home/ec2-user/efs/docker/NCBI/20200808/nr.dmnd"}

BASE_PATH=$(pwd)
LOCAL_QUERY_FASTA="${BASE_PATH}/${QUERY_PATH}"

LOCAL_QUERY_DIR="$(dirname ${LOCAL_QUERY_FASTA})"
QUERY_FILE_NAME="$(basename ${LOCAL_QUERY_FASTA})"

LOCAL_DB_DIR="$(dirname ${LOCAL_DB})"
DBNAME="$(basename ${LOCAL_DB})"
DB_PREFIX="$(basename ${DBNAME} .dmnd)"

# RESULTS_PATH=$(echo "${QUERY_PATH}" | sed -e "s/queries/results/")
RESULTS_PATH=${QUERY_PATH//"queries"/"results"}
LOCAL_RESULT_FILE="${BASE_PATH}/${RESULTS_PATH}"
LOCAL_RESULTS_DIR="$(dirname ${LOCAL_RESULT_FILE})"
RESULTS_PREFIX="$(basename ${LOCAL_RESULT_FILE} ${QUERY_EXT})_vs_${DB_PREFIX}"

RESULTS_FILE_NAME="${RESULTS_PREFIX}.dmnd_${TASK}.tsv"
TMP_DMND_OUTPUT="tmp__${RESULTS_PREFIX}"
DIAMOND_LOG="${RESULTS_PREFIX}.dmnd_${TASK}.log"

mkdir -p ${LOCAL_RESULTS_DIR}/${TMP_DMND_OUTPUT} ${LOCAL_RESULTS_DIR}/parallel_${TMP_DMND_OUTPUT}

DOCKER_DB_DIR="/app/db"
DOCKER_QUERIES_DIR="/app/queries"
DOCKER_RESULTS_DIR="/app/results"


docker container run --rm \
    -v ${LOCAL_DB_DIR}:${DOCKER_DB_DIR}:ro \
    -v ${LOCAL_QUERY_DIR}:${DOCKER_QUERIES_DIR}:ro \
    -v ${LOCAL_RESULTS_DIR}:${DOCKER_RESULTS_DIR}:rw \
    quay.io/biocontainers/diamond:2.0.2--h56fc30b_0 \
    diamond \
        ${TASK} \
            --verbose \
            --dbsize 1000000 \
            --threads ${THREADS} \
            --outfmt 6 \
            --tmpdir ${DOCKER_RESULTS_DIR}/${TMP_DMND_OUTPUT} \
            --parallel-tmpdir ${DOCKER_RESULTS_DIR}/parallel_${TMP_DMND_OUTPUT} \
            --db ${DOCKER_DB_DIR}/${DBNAME} \
            --out ${DOCKER_RESULTS_DIR}/${RESULTS_FILE_NAME} \
            --query ${DOCKER_QUERIES_DIR}/${QUERY_FILE_NAME} \
            --log &> ${LOCAL_RESULTS_DIR}/${DIAMOND_LOG}
