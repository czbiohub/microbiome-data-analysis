#!/bin/bash -x
# shellcheck disable=SC2086

set -e
set -u
set -o pipefail

DOCKER_IMAGE="ecogenomic/gtdbtk"
DOCKER_VERSION="1.1.1"

TASK="classify_wf"
THREADS=5
EXTENSION="fa"

INPUT_GENOMES_DIR="${1}"
OUTPUT_DIR="${2}"
LOCAL_DB_DIR="${3:-"/home/ec2-user/efs/docker/GTDB/SCv1/db/release89"}"
DOCKER_DB_PATH="/refdata"
DOCKER_DATA_PATH="/data"

docker container run --rm \
    --workdir "$(pwd)" \
    --volume "${LOCAL_DB_DIR}":"${DOCKER_DB_PATH}" \
    --volume "$(pwd)":"${DOCKER_DATA_PATH}" \
    ${DOCKER_IMAGE}:${DOCKER_VERSION} \
    ${TASK} \
        --cpus ${THREADS} \
        --extension ${EXTENSION} \
        --genome_dir ${DOCKER_DATA_PATH}/${INPUT_GENOMES_DIR} \
        --out_dir ${DOCKER_DATA_PATH}/${OUTPUT_DIR}