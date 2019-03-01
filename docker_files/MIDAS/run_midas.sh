#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154
# Requirements
#   DB (default): 36GB
#   memory: 2GB

set -e
set -u

LOCAL=/mnt
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
# database variable doesn't exist. 
#    exit

# get database from S3
# get path2db

COMMAND_PARAMS="-t ${coreNum} -d ${path2db} --remove-temp"

# Single end
if [[ -z ${fastq2} ]]; then
    COMMAND_PARAMS="${COMMAND_PARAMS} -1 ${fastq1}"
else
# Paired
    COMMAND_PARAMS="${COMMAND_PARAMS} -1 ${fastq1} -2 ${fastq2}"
fi

# Hard Trimming
if [[ -n $hardTrim ]]; then
    COMMAND_PARAMS="${COMMAND_PARAMS} --read_length ${hardTrim}"
fi

# Read subsets
if [[ -n ${subsetReads} ]]; then
    COMMAND_PARAMS="${COMMAND_PARAMS} -n ${subsetReads}"
fi

COMMAND="run_midas.py species ${OUTPUTDIR} ${COMMAND_PARAMS}"

eval "${COMMAND}"

aws s3 sync ${OUTPUTDIR} ${s3path}

rm -rf ${OUTPUTDIR}