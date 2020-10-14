#!/bin/bash -x
# shellcheck disable=SC2086 
# shellcheck disable=SC2046 
# USAGE: ls Fasta/ | parallel -j 10 --joblog hmmscan.joblog "bash run_hmmscan.sh {.} Proteins/{.}_protein.faa" &> run_hmmscan.log &

set -euo pipefail

DOCKER_IMAGE="quay.io/biocontainers/hmmer"
DOCKER_IMAGE_VERSION="3.3--he1b5a44_1"

PREFIX=${1}
FAA_FILE=${2}
NUM_CPU=6

DB_NAME="Pfam-A.hmm"
OUTPUT_DIR=${PREFIX}/nhmmscan
mkdir -p "${OUTPUT_DIR}"

docker container run --rm \
    --workdir $(pwd) \
    --volume $(pwd):$(pwd) \
    --volume /home/ec2-user/efs/docker/PFam/v33.1/PfamA:$(pwd)/db \
    ${DOCKER_IMAGE}:${DOCKER_IMAGE_VERSION} \
    nhmmscan \
        -o ${OUTPUT_DIR}/${PREFIX}.out.txt \
        --cpu ${NUM_CPU} \
        --max \
        --tblout ${OUTPUT_DIR}/${PREFIX}.tblout.txt \
        --domtblout ${OUTPUT_DIR}/${PREFIX}.domtblout.txt \
        --pfamtblout ${OUTPUT_DIR}/${PREFIX}.pfamtblout.txt \
        --acc \
        --notextw \
        --cut_tc \
        db/${DB_NAME} ${FAA_FILE} &> ${OUTPUT_DIR}/hmmscan.log

aws s3 sync --quiet ${OUTPUT_DIR} s3://czbiohub-microbiome/Sunit_Jain/Kanika/hydrogenotrophy/PFam_Search/${OUTPUT_DIR}
