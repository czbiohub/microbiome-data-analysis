#!/bin/bash -x
# shellcheck disable=SC2034
# shellcheck disable=SC2086
# shellcheck disable=SC2016

set -e
set -u
set -o pipefail

logger() {
  echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')]: $*" >&2
}

DOCKER_IMAGE=quay.io/biocontainers/magpurify
DOCKER_VERSION=2.1.2--py38h5ca1d4c_0

THREADS=50

INPUT_GENOME="${1}"
BAM_DIR="${2%/}"
OUTPUT_DIR="${3%/}"
CLEANED_OUTPUT="${OUTPUT_DIR}/$(basename ${1%.*}).mag_purify.fasta"

if [ -d ${OUTPUT_DIR} ]; then
    logger "Output directory already exists. Either delete it or choose another."
    exit 1
fi

# mash sketch -l example/ref_genomes.list -o example/ref_genomes
REF_GENOMES_MASH=${4:-"/home/ec2-user/efs/docker/UHGGdb/20200226_isolates/uhggdb_isolate_genomes.msh"} # ref_genomes.msh
REF_GENOMES_MASH=${REF_GENOMES_MASH%/}
LOCAL_DB_DIR="${5:-"/home/ec2-user/efs/docker/MAGpurify/MAGpurify-db-v1.0"}"
LOCAL_DB_DIR=${LOCAL_DB_DIR%/}

DOCKER_DB_PATH="$(pwd)/db"
DOCKER_MASH_PATH="$(pwd)/mash"

REF_GENOMES_MASH_FILE=$(basename ${REF_GENOMES_MASH})
REF_GENOMES_MASH_DIR=$(dirname ${REF_GENOMES_MASH})

# NOTE: the following variables are VERY SPECIFIC to how the mash sketch was created.
# Since the input to mash had a certain path, MAGpurify looks for the genome at the
# same path in the container. Thus another volume needs to be created.
# NOTE2: Cannot create mash of gz file, since the file is passed to blastn as is.
DOCKER_REF_PATH="$(pwd)/fasta"
REF_GENOMES_FASTA_DIR=${REF_GENOMES_MASH_DIR}/fasta

DOCKER_CMD_PREFIX='docker container run --rm --workdir "$(pwd)" --volume ${REF_GENOMES_FASTA_DIR}:${DOCKER_REF_PATH} --volume ${LOCAL_DB_DIR}:${DOCKER_DB_PATH} --volume ${REF_GENOMES_MASH_DIR}:${DOCKER_MASH_PATH} --volume "$(pwd)":"$(pwd)" --env MAGPURIFYDB=${DOCKER_DB_PATH} ${DOCKER_IMAGE}:${DOCKER_VERSION} magpurify'

# These modules support multi-threading
for TASK in phylo-markers clade-markers known-contam; do
    # magpurify phylo-markers example/test.fna example/output
    # magpurify clade-markers example/test.fna example/output
    # magpurify known-contam example/test.fna example/output
    eval ${DOCKER_CMD_PREFIX} \
        ${TASK} \
            --threads ${THREADS} \
            ${INPUT_GENOME} \
            ${OUTPUT_DIR} || \
    logger "An error occurred while trying to run \'${TASK}\' analysis. Moving on ..."
done

# These modules do not support multi-threading
for TASK in tetra-freq gc-content; do
    # magpurify tetra-freq example/test.fna example/output
    # magpurify gc-content example/test.fna example/output
    eval ${DOCKER_CMD_PREFIX} \
        ${TASK} \
            ${INPUT_GENOME} \
            ${OUTPUT_DIR} || \
    logger "An error occurred while trying to run \'${TASK}\' analysis. Moving on ..."
done

# magpurify conspecific example/test.fna example/output example/ref_genomes.msh
eval ${DOCKER_CMD_PREFIX} \
    conspecific \
        --threads ${THREADS} \
        ${INPUT_GENOME} \
        ${OUTPUT_DIR} \
        ${DOCKER_MASH_PATH}/${REF_GENOMES_MASH_FILE} || \
logger "An error occurred while trying to run 'conspecific' analysis. Moving on ..."

# magpurify coverage example/test.fna example/output BAM/sample_1.bam BAM/sample_2.bam BAM/sample_3.bam
eval ${DOCKER_CMD_PREFIX} \
    coverage \
        --threads ${THREADS} \
        ${INPUT_GENOME} \
        ${OUTPUT_DIR} \
        ${BAM_DIR}/*.bam || \
logger "An error occurred while trying to run 'coverage' analysis. Moving on ..."

# magpurify clean-bin example/test.fna example/output example/test_cleaned.fna
eval ${DOCKER_CMD_PREFIX} \
    clean-bin \
        ${INPUT_GENOME} \
        ${OUTPUT_DIR} \
        ${CLEANED_OUTPUT}
logger "Done!"