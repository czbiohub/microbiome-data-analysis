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

CHECKM_DOCKER_IMAGE="quay.io/biocontainers/checkm-genome"
CHECKM_DOCKER_VERSION="1.1.2--py_1"

GTDB_DOCKER_IMAGE="ecogenomic/gtdbtk"
GTDB_DOCKER_VERSION="1.2.0"

BBMAP_DOCKER_IMAGE="quay.io/biocontainers/bbmap"
BBMAP_DOCKER_VERSION="38.79--h516909a_0"

BARRNAP_DOCKER_IMAGE="quay.io/biocontainers/barrnap"
BARRNAP_DOCKER_VERSION="0.9--3"

QUAST_DOCKER_IMAGE="quay.io/biocontainers/quast"
QUAST_DOCKER_VERSION="5.0.2--1"

###############################################################################

SAMPLE_NAME=${1%/}
BIN_DIR=${2%/}
BIN_FASTA_EXT="${3:-"fa"}"
LOCAL_GTDB_DIR="${4:-"/home/ec2-user/efs/docker/GTDB/SCv1/db/release89"}"

STATS_DIR="${SAMPLE_NAME}/stats"
CHECKM_OUTPUT_DIR="${SAMPLE_NAME}/checkm"
CHECKM_OUTPUT_FILE="${CHECKM_OUTPUT_DIR}/${SAMPLE_NAME}.tsv"

GTDB_OUTPUT_DIR="${SAMPLE_NAME}/GTDBtk"
GTDB_DB_PATH="/refdata"
GTDB_DATA_PATH="/data"

BARRNAP_OUTPUT_DIR="${SAMPLE_NAME}/barrnap"
QUAST_OUTPUT_DIR="${SAMPLE_NAME}/quast"

mkdir -p ${STATS_DIR} ${CHECKM_OUTPUT_DIR} ${GTDB_OUTPUT_DIR} ${BARRNAP_OUTPUT_DIR} ${QUAST_OUTPUT_DIR}

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
        ${BIN_DIR}/*.${BIN_FASTA_EXT} 1> ${STATS_DIR}/fasta_stats.txt
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
            ${BIN_DIR} \
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
        --genome_dir ${GTDB_DATA_PATH}/${BIN_DIR} \
        --out_dir ${GTDB_DATA_PATH}/${GTDB_OUTPUT_DIR}

###############################################################################

# Run quast on each bin to predict the rRNA number 
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd)":"$(pwd)" \
    ${QUAST_DOCKER_IMAGE}:${QUAST_DOCKER_VERSION} \
	quast.py \
	-t ${THREADS} \
	--rna-finding \
	-o ${QUAST_OUTPUT_DIR} \
	${BIN_DIR}/*.${BIN_FASTA_EXT} 

# # Keep Genes and Proteins
# mkdir -p ${SAMPLE_NAME}/orfs/{genes,proteins}
# find ${GTDB_OUTPUT_DIR} -name '*_protein.fna' | parallel "cp {} ${SAMPLE_NAME}/orfs/genes/"
# find ${GTDB_OUTPUT_DIR} -name '*_protein.fna' | parallel "cp {} ${SAMPLE_NAME}/orfs/proteins/"

# sudo chmod a+wx ${GTDB_OUTPUT_DIR}/*/intermediate_results
# sudo rm -rf ${GTDB_OUTPUT_DIR}/*/intermediate_results

# grep -c '>' ${SAMPLE_NAME}/orfs/genes/*_protein.fna | \
#     sed -e 's/_protein.fna:/\t/' -e "s#${SAMPLE_NAME}/orfs/genes/##" > ${SAMPLE_NAME}/orfs/num_genes.txt
###############################################################################

# Compare sketch to NCBI
find ${BIN_DIR} -name "${SAMPLE_NAME}.*.${BIN_FASTA_EXT}" |\
parallel -j ${THREADS} \
    "docker container run --rm \
        --workdir $(pwd) \
        --volume $(pwd):$(pwd) \
        ${BBMAP_DOCKER_IMAGE}:${BBMAP_DOCKER_VERSION} \
        sendsketch.sh \
            in={} \
            out={.}.sketch \
            overwrite=t"

##############################################################################################################################################################

# Extract rrna
find ${BIN_DIR} -name "${SAMPLE_NAME}.*.${BIN_FASTA_EXT}" |\
parallel -j ${THREADS} \
    "docker container run --rm \
        --workdir $(pwd) \
        --volume $(pwd):$(pwd) \
        ${BARRNAP_DOCKER_IMAGE}:${BARRNAP_DOCKER_VERSION} \
        barrnap --threads 1 \
            -o ${BARRNAP_OUTPUT_DIR}/{/.}.rrna.${BIN_FASTA_EXT} {} > ${BARRNAP_OUTPUT_DIR}/{/.}.rrna.gff"