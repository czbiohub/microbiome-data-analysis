#!/bin/bash -x
set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

# S3INPUTPATH="s3://czb-seqbot/fastqs/190721_NB501938_0142_AHFFCYBGXB/Com2_0p1_Btmutpool_rep1_T12_S10_R1_001.fastq.gz"
# S3OUTPUTPATH="s3://czbiohub-microbiome/Sunit_Jain/BarSeq/FitnessBrowser/Com2_0p1_Btmutpool_rep1_T12_TEST"

echo "${PATH}"
LOCAL=$(pwd)
DEFAULT_CODE_DIR="/mnt"
SCRIPT_DIR="${DEFAULT_CODE_DIR}/FitnessBrowser/bin"
DEFAULT_INDEX_DIR="${DEFAULT_CODE_DIR}/FitnessBrowser/primers"
DEFAULT_INDEX_FILE="barseq3.index2"

S3INPUTPATH=${S3INPUTPATH%/}
S3OUTPUTPATH=${S3OUTPUTPATH%/}
MIN_QUALITY=${MIN_QUALITY:-20}
INDEX_DIR=${INDEX_DIR:-$DEFAULT_INDEX_DIR}
INDEX_FILE=${INDEX_FILE:-$DEFAULT_INDEX_FILE}

INDEX_FILE_LOC="${INDEX_DIR}/${INDEX_FILE}"

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
FIBO_OUTPUT="${LOCAL_OUTPUT}/fibo/codes"
LOCAL_FASTQ=$(basename "${S3INPUTPATH}")
PREFIX=$(basename "${S3OUTPUTPATH}")

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${FIBO_OUTPUT}"

trap '{aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}";
    rm -rf ${OUTPUTDIR} ; 
    exit 255; }' 1 

aws s3 cp --quiet "${S3INPUTPATH}" "${OUTPUTDIR}/${LOCAL_FASTQ}"

ls "${LOCAL}"
du -sh "${LOCAL}"

zcat "${OUTPUTDIR}/${LOCAL_FASTQ}" |\
    perl "${SCRIPT_DIR}/MultiCodes.pl" \
        -primers "${INDEX_FILE_LOC}" \
        -minQuality "${MIN_QUALITY}" \
        -out "${FIBO_OUTPUT}/${PREFIX}" |\
        tee -a "${LOG_DIR}/${PREFIX}.log"

ls "${LOCAL}"
du -sh "${LOCAL}"
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > "${LOCAL_OUTPUT}/job.complete"
echo "Live long and prosper" >> "${LOCAL_OUTPUT}/job.complete"
############################ PEACE! ################################
## Sync output
aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
