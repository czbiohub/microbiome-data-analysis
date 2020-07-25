#!/bin/bash -x

set -e
set -u
set -o pipefail

LOG_FILE="${1}"
TIMESTAMP=$( date +"%Y%m%d_%H%M%S" )
STATUS_FILE="${TIMESTAMP}_${LOG_FILE%.*}.status.csv"

echo "jobName,jobId,status,logStreamName" > ${STATUS_FILE}
cat ${LOG_FILE} | \
    jq -r .jobId | \
    parallel "aegea batch describe {} | jq -r '[.jobName,.jobId,.status,.container.logStreamName] | @csv'" >> \
    ${STATUS_FILE}