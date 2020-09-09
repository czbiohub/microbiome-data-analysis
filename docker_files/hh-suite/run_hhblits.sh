#!/bin/bash -x
# USAGE: ls Fasta/ | parallel -j 10 --joblog hhblits.joblog "bash run_hhblits.sh {.} Proteins/{.}_protein.faa" &> run_hhblits.log &
set -euo pipefail

DOCKER_IMAGE="soedinglab/hh-suite"
DOCKER_IMAGE_VERSION="latest"

PREFIX=${1}
FAA_FILE=${2}

NUM_JOBS=4
PROC_PER_JOB=4

TEMP_SEQ_DIR="tmp__${PREFIX}"
OUTPUT_DIR=${PREFIX}/hhblits
mkdir -p "${TEMP_SEQ_DIR}" "${OUTPUT_DIR}"

FAA_COPY=${PREFIX}.faa

cp ${FAA_FILE} ${TEMP_SEQ_DIR}/${FAA_COPY}

# Given a fasta file split it into N files, where N == number of sequences.
docker container run --rm \
    --workdir $(pwd)/${TEMP_SEQ_DIR} \
    --volume $(pwd)/${TEMP_SEQ_DIR}:$(pwd)/${TEMP_SEQ_DIR} \
    soedinglab/hh-suite \
        splitfasta.pl ${FAA_COPY}

rm ${TEMP_SEQ_DIR}/${FAA_COPY}

# Find sequences in UniRef30 similar to our sequence of interest and create an MSA
find ${TEMP_SEQ_DIR} -name "*.seq" \
    | parallel -j ${NUM_JOBS} --joblog ${PREFIX}.msa.joblog \
    "docker container run --rm \
        --workdir $(pwd) \
        --volume $(pwd)/:$(pwd) \
        --volume $(pwd)/uniref:/db \
        ${DOCKER_IMAGE}:${DOCKER_IMAGE_VERSION} \
        hhblits \
            -i {} \
            -d /db/UniRef30_2020_03 \
            -oa3m ${TEMP_SEQ_DIR}/{/.}.a3m  \
            -cpu ${PROC_PER_JOB} \
            -n 2" &> ${OUTPUT_DIR}/hhblits.UniRef30_2020_03.log

# Use the MSA to find similar profiles in PFamA
find ${TEMP_SEQ_DIR} -name "*.a3m" \
    | parallel -j ${NUM_JOBS} --joblog ${PREFIX}.joblog \
    "docker container run --rm \
        --workdir $(pwd) \
        --volume $(pwd)/:$(pwd) \
        --volume $(pwd)/PfamA_v32:/db \
        ${DOCKER_IMAGE}:${DOCKER_IMAGE_VERSION} \
        hhblits \
            -cpu ${PROC_PER_JOB} \
            -n 2 \
            -i {} \
            -d /db/pfam \
            -o ${OUTPUT_DIR}/{/.}.hhr \
            -blasttab ${OUTPUT_DIR}/{/.}.blasttab.tsv" &> ${OUTPUT_DIR}/hhblits.pfamA.log

rm -rf ${TEMP_SEQ_DIR}

aws s3 sync --quiet \
    ${OUTPUT_DIR} \
    s3://czbiohub-microbiome/Sunit_Jain/Kanika/hydrogenotrophy/PFam_Search/${OUTPUT_DIR}
