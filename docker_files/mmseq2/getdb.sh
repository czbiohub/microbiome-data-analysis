#!/bin/bash -x
set -euo pipefail

# For valid names, refer: 
# https://github.com/soedinglab/MMseqs2/wiki#downloading-databases
DB_NAME=${1} 

SAFE_DB_NAME=$(echo "${DB_NAME}" | tr "/" "_")
SAFE_DB_NAME=${SAFE_DB_NAME,,}

mkdir -p "${SAFE_DB_NAME}"

docker container run --rm \
    --volume "$(pwd)/${SAFE_DB_NAME}":/app \
    soedinglab/mmseqs2 \
    mmseqs \
        databases \
            ${DB_NAME} \
            /app/${SAFE_DB_NAME} \
            /app/${SAFE_DB_NAME}_tmp