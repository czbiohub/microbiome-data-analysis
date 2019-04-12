#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-16};
LOCAL_DB_PATH=${LOCAL}/databases
DB_VERSION=${DB_VERSION:-"v1"}
S3DBPATH=s3://czbiohub-microbiome/ReferenceDBs/kraken2/minikraken2/${DB_VERSION}
# fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/Bracken_Test/Dorea-longicatena-DSM-13814

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
BBDUK_OUTPUT="${LOCAL_OUTPUT}/BBDuk"
KRAKEN2_OUTPUT="${LOCAL_OUTPUT}/kraken2"
BRACKEN_OUTPUT="${LOCAL_OUTPUT}/bracken"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

# Create intermediate directories
mkdir -p "${OUTPUTDIR}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${LOCAL_DB_PATH}"
# Create output directories that will be syncked with S3
mkdir -p "${LOCAL_OUTPUT}" "${BBDUK_OUTPUT}" "${KRAKEN2_OUTPUT}" "${BRACKEN_OUTPUT}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1 

hash_kmer=${hash_kmer:-51}

cd "${LOCAL}" || exit

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

###############################################################################

# Constant definitions for bbduk
adapterFile="adapters,phix"
trimQuality=${trimQuality:-25}
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

# Use bbduk to trim reads, -eoom exits when out of memory
bbduk.sh -Xmx16g tbo -eoom hdist=1 qtrim=rl ktrim=r \
    in1="${RAW_FASTQ}/read1.fastq.gz" \
    in2="${RAW_FASTQ}/read2.fastq.gz" \
    out1="${QC_FASTQ}/read1_trimmed.fastq.gz" \
    out2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
    ref=${adapterFile} \
    k="${kmer_value}" \
    mink="${min_kmer_value}" \
    trimq="${trimQuality}" \
    minlen="${minLength}" \
    refstats="${BBDUK_OUTPUT}/adapter_trimming_stats_per_ref.txt" |\
    tee -a ${LOG_DIR}/bbduk.log
###############################################################################

# Get database
aws s3 sync --quiet ${S3DBPATH} ${LOCAL_DB_PATH}
# tar xvzf --directory ${LOCAL_DB_PATH}/ ${LOCAL_DB_PATH}/minikraken2_v1_8GB.tgz
# KRAKEN2_DB=${LOCAL_DB_PATH}/minikraken2_${DB_VERSION}_8GB

# Run Kraken2
kraken2 --gzip-compressed --paired\
    --db=${LOCAL_DB_PATH} \
    --threads ${coreNum} \
    --report ${KRAKEN2_OUTPUT}/${SAMPLE_NAME}.kreport2.tsv \
    --output ${OUTPUTDIR}/${SAMPLE_NAME}.kraken2 \
    "${QC_FASTQ}/read1_trimmed.fastq.gz" \
    "${QC_FASTQ}/read2_trimmed.fastq.gz" |\
    tee -a ${LOG_DIR}/kraken2.log
    
###############################################################################

READ_LEN=${READ_LEN:-150}
CLASSIFICATION_LEVEL=${CLASSIFICATION_LEVEL:-"S"}
THRESHOLD=${THRESHOLD:-100}

# Run Bracken
bracken \
    -d ${LOCAL_DB_PATH} \
    -i ${KRAKEN2_OUTPUT}/${SAMPLE_NAME}.kreport2.tsv \
    -o ${BRACKEN_OUTPUT}/${SAMPLE_NAME}.bracken.tsv \
    -r ${READ_LEN} \
    -l ${CLASSIFICATION_LEVEL} \
    -t ${THRESHOLD} |\
    tee -a ${LOG_DIR}/kraken2.log

# Alternative:
# python estimate_abundance.py \
#     -i ${SAMPLE_NAME}.kreport \
#     -k ${KRAKEN_DB}/database${READ_LEN}mers.kmer_distrib \
#     -o ${SAMPLE_NAME}.bracken \
#     -l ${CLASSIFICATION_LEVEL} \
#     -t ${THRESHOLD}
###############################################################################

aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
