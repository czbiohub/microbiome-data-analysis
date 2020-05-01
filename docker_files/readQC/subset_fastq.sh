#!/bin/bash -x

set -e
set -u
set -o pipefail

FASTQ=${1}
OUTPUT_PREFIX=${2}
NUM_SEQS=${3}
subset_fastq() {
    zcat "${1}" | \
        paste - - - - | \
        shuf -n ${3} | \
        awk '{printf "%s\n%s\n", $1, $2}' | \
        sed 's/@/>/' > "${2}_shuffled.n${3}.fasta"
}

subset_fastq "${FASTQ}" "${OUTPUT_PREFIX}" "${NUM_SEQS}"
