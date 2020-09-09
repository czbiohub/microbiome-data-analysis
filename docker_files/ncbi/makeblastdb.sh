#!/bin/bash -x

set -euo pipefail

DBFASTA="${1}"
DBNAME="${2}"
DBTYPE="${3:-'nucl'}"

mkdir -p blastdb_custom fasta
mv ${DBFASTA} fasta/
DB_FASTA_NAME="$(basename ${DBFASTA})"

docker run --rm \
    -v $(pwd)/blastdb_custom:/blast/blastdb_custom:rw \
    -v $(pwd)/fasta:/blast/fasta:ro \
    -w /blast/blastdb_custom \
    ncbi/blast \
        makeblastdb \
            -in /blast/fasta/${DB_FASTA_NAME} \
            -dbtype ${DBTYPE} \
            -out ${DBNAME}