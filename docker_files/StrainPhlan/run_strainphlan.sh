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

# S3DBPATH=s3://czbiohub-microbiome/ReferenceDBs/Sourmash/Genbank/2018-03/genbank-d2-k51.tar.gz
# fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/Sourmash_Test/Dorea-longicatena-DSM-13814

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
INPUT_FASTQ="${OUTPUTDIR}/raw_fastq"
CAT_FASTQ="${OUTPUTDIR}/cat_fastq"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
MASH_OUTPUT="${LOCAL_OUTPUT}/mash"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${CAT_FASTQ}" "${INPUT_FASTQ}" "${MASH_OUTPUT}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1 

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${INPUT_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${INPUT_FASTQ}/read2.fastq.gz"
cat ${INPUT_FASTQ}/read1.fastq.gz ${INPUT_FASTQ}/read2.fastq.gz > ${CAT_FASTQ}/concat.fastq.gz

# Constant definitions
num_sketches=${num_sketches:-1000000}
hash_kmer=${hash_kmer:-32}
min_kmer_value=${min_kmer_value:-5}

# Compute the signature of the input fastq file
mash sketch \
    -r \
    -p ${coreNum} \
    -k ${hash_kmer} \
    -s ${num_sketches}\
    -m ${min_kmer_value} \
    -o ${MASH_OUTPUT}/${SAMPLE_NAME}.msh \
    ${CAT_FASTQ}/concat.fastq.gz

aws s3 sync --quiet "${MASH_OUTPUT}" "${S3OUTPUTPATH}"
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${MASH_OUTPUT}/job.complete
############################ PEACE! ################################
aws s3 sync --quiet "${MASH_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
