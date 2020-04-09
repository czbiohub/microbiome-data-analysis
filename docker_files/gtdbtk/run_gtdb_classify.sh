#!/bin/bash -x
# shellcheck disable=SC2086

set -e
set -u
set -o pipefail

DOCKER_IMAGE=sunitjain/gtdbtk
DOCKER_VERSION=20200107150123

# DOCKER_IMAGE=quay.io/biocontainers/gtdbtk
# DOCKER_VERSION=1.0.2--py_3

TASK="classify_wf"
THREADS=50
EXTENSION="fa"

INPUT_GENOMES_DIR="${1}"
OUTPUT_DIR="${2}"
LOCAL_DB_DIR="${3:-"/home/ec2-user/efs/docker/GTDB/SCv1/db/release89"}"
DOCKER_DB_PATH="/opt/conda/share/gtdbtk-0.3.2/db/"

docker container run \
   --workdir "$(pwd)" \
   --volume ${LOCAL_DB_DIR}:${DOCKER_DB_PATH} \
   --volume "$(pwd)":"$(pwd)" \
   --env GTDBTK_DATA_PATH=${DOCKER_DB_PATH} \
   ${DOCKER_IMAGE}:${DOCKER_VERSION} \
   gtdbtk \
        ${TASK} \
            --cpus ${THREADS} \
            --extension ${EXTENSION} \
            --genome_dir ${INPUT_GENOMES_DIR} \
            --out_dir ${OUTPUT_DIR}