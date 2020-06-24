#!/usr/bin/bash -x
# shellcheck disable=SC2086

set -e
set -u
set -o pipefail

THREADS=15
# PP_THREADS=$(("${THREADS}/2"))
PP_THREADS=7

SEQKIT_DOCKER_IMAGE="quay.io/biocontainers/seqkit"
SEQKIT_DOCKER_VERSION="0.12.0--0"

METABAT_DOCKER_IMAGE="metabat/metabat"
METABAT_DOCKER_VERSION="latest"

CHECKM_DOCKER_IMAGE="quay.io/biocontainers/checkm-genome"
CHECKM_DOCKER_VERSION="1.1.2--py_1"

GTDB_DOCKER_IMAGE="ecogenomic/gtdbtk"
GTDB_DOCKER_VERSION="1.2.0"

BBMAP_DOCKER_IMAGE="quay.io/biocontainers/bbmap"
BBMAP_DOCKER_VERSION="38.79--h516909a_0"

BARRNAP_DOCKER_IMAGE="quay.io/biocontainers/barrnap"
BARRNAP_DOCKER_VERSION="0.9--3"

###############################################################################

SAMPLE_NAME=${1}
ASSEMBLED_FASTA=${2}
BAM_DIR=${3%/}

BIN_FASTA_EXT="${4:-"fa"}"
LOCAL_GTDB_DIR="${5:-"/home/ec2-user/efs/docker/GTDB/SCv1/db/release89"}"

STATS_DIR="${SAMPLE_NAME}/stats"
METABAT_OUTPUT_DIR="${SAMPLE_NAME}/metabat"
METABAT_DEPTH_FILE="${METABAT_OUTPUT_DIR}/${SAMPLE_NAME}.depth.txt"
mkdir -p ${METABAT_OUTPUT_DIR}

CHECKM_OUTPUT_DIR="${SAMPLE_NAME}/checkm"
CHECKM_OUTPUT_FILE="${CHECKM_OUTPUT_DIR}/${SAMPLE_NAME}.tsv"

GTDB_OUTPUT_DIR="${SAMPLE_NAME}/GTDBtk"
GTDB_DB_PATH="/refdata"
GTDB_DATA_PATH="/data"

BARRNAP_OUTPUT_DIR="${SAMPLE_NAME}/barrnap"

mkdir -p ${STATS_DIR} ${METABAT_OUTPUT_DIR} ${CHECKM_OUTPUT_DIR} ${GTDB_OUTPUT_DIR} ${BARRNAP_OUTPUT_DIR}
###############################################################################
# Run SeqKit for BAM stats
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd)":"$(pwd)" \
    ${SEQKIT_DOCKER_IMAGE}:${SEQKIT_DOCKER_VERSION} \
    seqkit \
        bam \
        -s \
        -j ${THREADS} \
        ${BAM_DIR}/*bam 1> ${STATS_DIR}/bam_stats.txt

###############################################################################
# Run Metabat2
## jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd)":"$(pwd)" \
    ${METABAT_DOCKER_IMAGE}:${METABAT_DOCKER_VERSION} \
    jgi_summarize_bam_contig_depths \
        --outputDepth ${METABAT_DEPTH_FILE} \
        ${BAM_DIR}/*bam

## metabat2 -i assembly.fasta -a depth.txt -o bins_dir/bin
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd)":"$(pwd)" \
    ${METABAT_DOCKER_IMAGE}:${METABAT_DOCKER_VERSION} \
    metabat2 \
        --seed 1712 \
        -t ${THREADS} \
        -i ${ASSEMBLED_FASTA} \
        -a ${METABAT_DEPTH_FILE} \
        -o ${METABAT_OUTPUT_DIR}/${SAMPLE_NAME}

###############################################################################
# Run SeqKit for fasta stats
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd)":"$(pwd)" \
    ${SEQKIT_DOCKER_IMAGE}:${SEQKIT_DOCKER_VERSION} \
    seqkit \
        stats \
        -T \
        -j ${THREADS} \
        ${ASSEMBLED_FASTA} ${METABAT_OUTPUT_DIR}/*.${BIN_FASTA_EXT} 1> ${STATS_DIR}/fasta_stats.txt

###############################################################################
# Run CheckM on each bin from metaBAT
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd)":"$(pwd)" \
    ${CHECKM_DOCKER_IMAGE}:${CHECKM_DOCKER_VERSION} \
    checkm \
        lineage_wf \
            --nt \
            --ali \
            --tab_table \
            -t ${THREADS} \
            -x ${BIN_FASTA_EXT} \
            --pplacer_threads ${PP_THREADS} \
            -f ${CHECKM_OUTPUT_FILE} \
            ${METABAT_OUTPUT_DIR} \
            ${CHECKM_OUTPUT_DIR}

###############################################################################

# place each bin on a tree (GTDB Classify).
# --env GTDBTK_DATA_PATH=${DOCKER_DB_PATH} \
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "${LOCAL_GTDB_DIR}":"${GTDB_DB_PATH}" \
    --volume "$(pwd)":"${GTDB_DATA_PATH}" \
    ${GTDB_DOCKER_IMAGE}:${GTDB_DOCKER_VERSION} \
    classify_wf \
        --cpus ${THREADS} \
        --extension ${BIN_FASTA_EXT} \
        --genome_dir ${GTDB_DATA_PATH}/${METABAT_OUTPUT_DIR} \
        --out_dir ${GTDB_DATA_PATH}/${GTDB_OUTPUT_DIR}

###############################################################################

# Compare sketch to NCBI
find ${METABAT_OUTPUT_DIR} -name "${SAMPLE_NAME}.*.${BIN_FASTA_EXT}" |\
parallel -j ${THREADS} \
    "docker container run --rm \
        --workdir $(pwd) \
        --volume $(pwd):$(pwd) \
        ${BBMAP_DOCKER_IMAGE}:${BBMAP_DOCKER_VERSION} \
        sendsketch.sh \
            in={} \
            out={.}.sketch \
            overwrite=t"

###############################################################################

# Extract rrna
find ${METABAT_OUTPUT_DIR} -name "${SAMPLE_NAME}.*.${BIN_FASTA_EXT}" |\
parallel -j ${THREADS} \
    "docker container run --rm \
        --workdir $(pwd) \
        --volume $(pwd):$(pwd) \
        ${BARRNAP_DOCKER_IMAGE}:${BARRNAP_DOCKER_VERSION} \
        barrnap --threads 1 \
            -o ${BARRNAP_OUTPUT_DIR}/{/.}.rrna.${BIN_FASTA_EXT} {} > ${BARRNAP_OUTPUT_DIR}/{/.}.rrna.gff"
