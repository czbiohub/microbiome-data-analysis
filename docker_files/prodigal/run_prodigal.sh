#!/bin/bash -x

set -euo pipefail

NAME=${1}
S3CONTIG=${2}
S3OUTPUTPATH=${3:-"s3://czbiohub-microbiome/Sonnenburg_Lab/InfantMicrobiome/Prodigal"}
S3OUTPUTPATH=${S3OUTPUTPATH%/}

# Setup
PRODIGAL_DOCKER_IMAGE="quay.io/biocontainers/prodigal"
PRODIGAL_IMAGE_VERSION="2.6.3--h516909a_2"
RENAME_SCRIPT="$(pwd)/scripts/remove_partials.py"

WORKDIR="$(pwd)/${NAME}"
INDIR="${WORKDIR}/IN"
OUTDIR="${WORKDIR}/OUT"

# Local fastq files
LOCAL_CONTIG=${INDIR}/contigs.fasta

# Download files from S3
aws s3 cp ${S3CONTIG} ${LOCAL_CONTIG} --quiet

# Create output directory
mkdir -p ${OUTDIR}

# Run Prodigal
docker container run --rm \
    --workdir "${WORKDIR}" \
    --volume "${INDIR}:${INDIR}" \
    --volume "${OUTDIR}:${OUTDIR}" \
    ${PRODIGAL_DOCKER_IMAGE}:${PRODIGAL_IMAGE_VERSION} \
    prodigal \
        -p meta \
        -f gff \
        -a ${OUTDIR}/proteins.faa \
        -d ${OUTDIR}/genes.fna \
        -o ${OUTDIR}/genes.gff \
        -s ${OUTDIR}/gene_scores.txt \
        -i ${INDIR}/contigs.fasta &> ${OUTDIR}/${NAME}.prodigal.log


python ${RENAME_SCRIPT} ${OUTDIR}/genes.fna ${NAME} ${OUTDIR} "partials_removed.genes.fna" &> ${OUTDIR}/${NAME}.fna_partials_removed.log
python ${RENAME_SCRIPT} ${OUTDIR}/proteins.faa ${NAME} ${OUTDIR} "partials_removed.proteins.faa" &> ${OUTDIR}/${NAME}.faa_partials_removed.log

# Remove original files
rm -rf ${INDIR}

# Sync repaired files to S3
aws s3 sync ${OUTDIR} ${S3OUTPUTPATH}/${NAME} --quiet

# CleanUp local directory
rm -rf ${NAME}
