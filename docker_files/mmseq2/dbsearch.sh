#!/bin/bash -x

set -euo pipefail

QUERY="${1}"
DBNAME="${2}"
FULL_DB_DIR=${3:-"/home/ec2-user/efs/docker/MMSeq2"}

mkdir -p Query Result
cp ${QUERY} Query/

USE_DB="${DBNAME}/${DBNAME}"
if [[ -f ${FULL_DB_DIR}/${DBNAME}/${DBNAME}.index ]]; then
    echo "DB already exists. Yay!"
else
    cp ${DBNAME} ${FULL_DB_DIR}/
    USE_DB="${DB_FILE}"
fi

QUERY_FILE=$(basename ${QUERY})
DB_FILE=$(basename ${DBNAME})
QUERY_BASE=$(basename ${QUERY_FILE} .fasta)
DB_BASE=$(basename ${DB_FILE} .fasta)

# Create Query DB
docker container run --rm \
    --volume "$(pwd)/Query":/app/input \
    --volume "$(pwd)/Result":/app/output \
    --volume "${FULL_DB_DIR}":/app/db \
    soedinglab/mmseqs2 \
    mmseqs \
        easy-search \
        /app/input/${QUERY_FILE} \
        /app/db/${USE_DB} \
        /app/output/${QUERY_BASE}_vs_${DB_BASE}.m8 \
        /app/tmp \
        -s 8 \
        --db-load-mode 3

# docker container run --rm \
#     --volume /home/ec2-user/efs/docker/Sonnenburg_Lab/InfantMicrobiome/MMSeqs2/Query:/app/input \
#     --volume /home/ec2-user/efs/docker/Sonnenburg_Lab/InfantMicrobiome/MMSeqs2/Result:/app/output \
#     --volume /home/ec2-user/efs/docker/MMSeq2:/app/db \
#     soedinglab/mmseqs2 \
#         mmseqs 
#             easy-search \
#             /app/input/infantMicrobiome_genes_rep_seq.fasta \
#             /app/db/nr/nr \
#             /app/output/infantMicrobiome_genes_rep_seq_vs_nr.m8 \
#             /app/tmp