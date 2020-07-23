#!/bin/bash -x
# NOTE: If running this script in parallel, make sure that each job is running
# in its own directory. RandomReads.sh creates an index directory that is always
# named ref. 
# If multiple jobs are run from same dir, they will overwrite this dir.

set -euo pipefail

REF_FASTA=${1}
COV=${2}
OUTPUTDIR=${3:-"${COV}x"}

MAX_SNPS=${MAX_SNPS:-0}
SNPS_ERR=${SNPS_ERR:-0}
SEQ_ERR=${SEQ_ERR:-0}

if [[ $COV -le 0 ]]; then
    exit 0
fi

mkdir -p "${OUTPUTDIR}"

SAMPLE_NAME=$(basename "${REF_FASTA}" .fna)
OUT_PREFIX=${OUTPUTDIR}/${SAMPLE_NAME}__covx${COV}_snp${MAX_SNPS}

BBTOOLS_IMAGE_NAME="quay.io/biocontainers/bbmap"
BBTOOLS_IMAGE_VERSION="38.79--h516909a_0"

docker container run --rm \
    --workdir $(pwd) \
    -v $(pwd):$(pwd) \
    ${BBTOOLS_IMAGE_NAME}:${BBTOOLS_IMAGE_VERSION} \
    randomreads.sh \
        maxsnps=${MAX_SNPS} \
        snprate=${SNPS_ERR} \
        adderrors=${SEQ_ERR} \
        paired=t \
        superflat=t \
        coverage=${COV} \
        minlength=140 \
        maxlength=140 \
        seed=1712 \
        out1=${OUT_PREFIX}.R1.fastq \
        out2=${OUT_PREFIX}.R2.fastq \
        ref=${REF_FASTA}

