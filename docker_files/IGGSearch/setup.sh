#!/bin/bash -x

set -e
set -u
set -o pipefail

LOCAL=$(pwd)
S3_IGGSEARCH_CODEBASE="s3://czbiohub-microbiome/ReferenceDBs/IGGdb/IGGsearch-master.zip"
IGGSEARCH_CODEBASE_DIR=$(basename ${S3_IGGSEARCH_CODEBASE} .zip)
IGG_CODE_PATH=${LOCAL}/src

aws s3 cp ${S3_IGGSEARCH_CODEBASE} "${IGG_CODE_PATH}"/

cd "${IGG_CODE_PATH}" && \
    unzip "${IGGSEARCH_CODEBASE_DIR}.zip" &&
    rm "${IGGSEARCH_CODEBASE_DIR}.zip" && \
    cd "${LOCAL}"

export PYTHONPATH="${IGG_CODE_PATH}/${IGGSEARCH_CODEBASE_DIR}/iggsearch"
export PATH="${IGG_CODE_PATH}/${IGGSEARCH_CODEBASE_DIR}:/opt/conda/bin:${PATH}"