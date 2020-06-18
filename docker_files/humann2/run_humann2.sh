#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-15};
LOCAL_DB_PATH=${LOCAL}/databases

# fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/Bowtie2_Test/Dorea-longicatena-DSM-13814
# S3DBPREFIX=s3://czbiohub-microbiome/Synthetic_Community/Genome_References/ncbi_fasta/Dorea-longicatena-DSM-13814-GCF_000154065.1_ASM15406v1

S3DBPATH=${S3DBPATH:-"s3://czbiohub-microbiome/ReferenceDBs/Humann2"}
DEFAULT_MP_DB_PATH="/opt/conda/bin/metaphlan_databases"

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
TMP_HM2_OUTPUT="${OUTPUTDIR}/tmp"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
HM2_OUTPUT="${LOCAL_OUTPUT}/Humann2"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})
FASTQ_FILES=${TMP_HM2_OUTPUT}/combined.fastq.gz

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}"
mkdir -p "${LOCAL_DB_PATH}" "${HM2_OUTPUT}" "${TMP_HM2_OUTPUT}"
trap '{ 
    aws s3 sync ${LOCAL_OUTPUT}/ ${S3OUTPUTPATH}/;
    rm -rf ${OUTPUTDIR} ; 
    exit 255; 
    }' 1 
# Copy QC-ed fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1_trimmed.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2_trimmed.fastq.gz"

# Concatenate
cat \
    ${RAW_FASTQ}/read1_trimmed.fastq.gz \
    ${RAW_FASTQ}/read2_trimmed.fastq.gz > \
    ${FASTQ_FILES}

## Get Databases ##
aws s3 sync ${S3DBPATH} ${LOCAL_DB_PATH} --quiet

## Set Databases ##
humann2_config --update database_folders utility_mapping ${LOCAL_DB_PATH}/utility_mapping
humann2_config --update database_folders nucleotide ${LOCAL_DB_PATH}/chocophlan
humann2_config --update database_folders protein ${LOCAL_DB_PATH}/uniref
METAPHLAN_BOWTIE_DB=${LOCAL_DB_PATH}/metaphlan2

# Set the default metaphlan2 database
rm -rf "${DEFAULT_MP_DB_PATH}" || exit 1
ln -s ${METAPHLAN_BOWTIE_DB} "${DEFAULT_MP_DB_PATH}"

JOB_LOG="${LOG_DIR}/${SAMPLE_NAME}.log"

## Humann2 the reads
humann2 \
    --verbose \
    --nucleotide-database ${LOCAL_DB_PATH}/chocophlan \
    --protein-database ${LOCAL_DB_PATH}/uniref \
    --o-log ${JOB_LOG} \
    --threads ${coreNum} \
    --output-basename ${SAMPLE_NAME} \
    --input ${FASTQ_FILES} \
    --output ${HM2_OUTPUT} \
    --memory-use maximum \
    --remove-temp-output

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
rm -rf "${OUTPUTDIR}"