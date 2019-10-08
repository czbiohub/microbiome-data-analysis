#!/bin/bash -x
set -e
set -u
set -o pipefail

START_TIME=$SECONDS
coreNum=${coreNum:-4}
export PATH="/opt/conda/bin:${PATH}"
export MC_CORES=${coreNum}

# S3INPUTPATH="s3://czb-seqbot/fastqs/190721_NB501938_0142_AHFFCYBGXB/Com2_0p1_Btmutpool_rep1_T12_S10_R1_001.fastq.gz"
# INDEX_NAME=""
# S3OUTPUTPATH="s3://czbiohub-microbiome/Sunit_Jain/BarSeq/FitnessBrowser/Com2_0p1_Btmutpool_rep1_T12_TEST"
# ORG_NAME

echo "${PATH}"
LOCAL=$(pwd)
DEFAULT_CODE_DIR="/mnt"
SCRIPT_DIR="${DEFAULT_CODE_DIR}/FitnessBrowser/bin"
DEFAULT_INDEX_DIR="${DEFAULT_CODE_DIR}/FitnessBrowser/primers"
DEFAULT_INDEX_FILE="barseq3.index2"
S3ASSETS="s3://czbiohub-microbiome/BarSeq/FitnessBrowser/Organisms"

S3INPUTPATH=${S3INPUTPATH%/}
S3OUTPUTPATH=${S3OUTPUTPATH%/}
MIN_QUALITY=${MIN_QUALITY:-20}
INDEX_DIR=${INDEX_DIR:-$DEFAULT_INDEX_DIR}
INDEX_FILE=${INDEX_FILE:-$DEFAULT_INDEX_FILE}

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
# FIBO_OUTPUT="${LOCAL_OUTPUT}/fibo/codes"
LOCAL_FASTQ=$(basename "${S3INPUTPATH}")
PREFIX=$(basename "${S3OUTPUTPATH}")

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${FIBO_OUTPUT}"

trap '{aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}";
    rm -rf ${OUTPUTDIR} ; 
    exit 255; }' EXIT
# Step 0: Download input fastqs and include indices intheir names.
aws s3 cp --quiet "${S3INPUTPATH}" "${OUTPUTDIR}/${LOCAL_FASTQ}"

# Step 0.1
# Rename to include ITXXX number.

ls "${LOCAL}"
du -sh "${LOCAL}"

# Step 1:
# Download some asset files and organize it in a directory like so:
mkdir -p "g/${ORG_NAME}"
aws s3 sync "${S3ASSETS}/${ORG_NAME}" "g/${ORG_NAME}"
# g
# └── ORG_NAME
#     ├── ORG_NAME.pool -> ../../assets/ORG_NAME.pool
#     ├── genes.GC -> ../../assets/ORG_NAME.GC
#     └── pool -> ORG_NAME.pool

# if genes.GC does not exist, create it using the following command.
# Requires a genomes file and a gff file.
# ${SCRIPT_DIR}/RegionGC.pl \
#     Price_Btheta_Genome.fna \
#     organism_Btheta_genes.txt > Price_Btheta_Genome.GC

INDICES="IT018:IT042:IT006:IT012:IT024:IT030:IT036:IT048:IT054:IT066:IT078:IT090"
INDEX_DESC="Time0:PBS:1:2:3:4:5:6:7:8:9:10"
# Step 2A
perl ${SCRIPT_DIR}/BarSeqTest.pl \
    -bs3 \
    -org "${ORG_NAME}" \
    -index "${INDICES}" \
    -desc "${INDEX_DESC}" \
    -fastqdir "${OUTPUTDIR}/${LOCAL_FASTQ}"

# OR
# Step 2B.1
# Do this for each sample in a set; maybe use parallel?
# zcat "${OUTPUTDIR}/${LOCAL_FASTQ}" |\
#     perl "${SCRIPT_DIR}/MultiCodes.pl" \
#         -bs3 \
#         -index "${INDEX_NAME}" \
#         -minQuality "${MIN_QUALITY}" \
#         -out "${FIBO_OUTPUT}/${PREFIX}" |\
#     tee -a "${LOG_DIR}/${PREFIX}.log"

# Step 2B.2
# ${SCRIPT_DIR}/BarSeqR.pl \
#     -org bTheta \
#     -exps g/bTheta/barseqtest/exps_table \
#     -pool g/bTheta/pool \
#     -indir g/bTheta/barseqtest \
#     -outdir g/bTheta/barseqtest \
#     -genes g/bTheta/genes.GC \
#     test

# Step 2B.3
# Rscript ${SCRIPT_DIR}/RunFEBA.R \
#     bTheta \
#     g/bTheta/barseqtest \
#     /mnt/FitnessBrowser/bin/.. > g/bTheta/barseqtest/log

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
