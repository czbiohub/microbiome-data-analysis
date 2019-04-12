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
METHOD=${METHOD:-"gather"}

S3DBPATH=s3://czbiohub-microbiome/ReferenceDBs/Sourmash/Genbank/2018-03/genbank-d2-k51.tar.gz
# fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/Sourmash_Test/Dorea-longicatena-DSM-13814

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
SOURMASH_OUTPUT="${LOCAL_OUTPUT}/sourmash"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${SOURMASH_OUTPUT}"
mkdir -p "${LOCAL_DB_PATH}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1 

hash_kmer=${hash_kmer:-51}
mkdir -p ${SOURMASH_OUTPUT}/compute ${SOURMASH_OUTPUT}/gather

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

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
    refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt" |\
    tee -a ${LOG_DIR}/bbduk.log

# Use khmer to remove low abundance reads that may be artifacts (optional step)
# trim-low-abund.py -C 3 -Z 18 -V -M 60e9 \
#     --summary-info tsv --quiet \
#     "${QC_FASTQ}/read1_trimmed.fastq.gz" \
#     "${QC_FASTQ}/read2_trimmed.fastq.gz"

# Compute the signature of the input fastq file
sourmash compute --scaled 1000 -k ${hash_kmer} \
    "${QC_FASTQ}/read1_trimmed.fastq.gz" \
    "${QC_FASTQ}/read2_trimmed.fastq.gz" \
    --track-abundance \
    --merge ${SAMPLE_NAME} \
    -o ${SOURMASH_OUTPUT}/compute/${SAMPLE_NAME}-reads.sig |\
    tee -a ${LOG_DIR}/sourmash_compute_reads_sig.log

# rm -rf read*_trimmed*.abundtrim
# Get database
LOCAL_DBNAME=$(basename ${S3DBPATH} .tar.gz)
SOURMASH_PARAMS=""

aws s3 cp --quiet ${S3DBPATH} ${LOCAL_DB_PATH}/
cd ${LOCAL_DB_PATH} || exit
tar xzf ${LOCAL_DBNAME}.tar.gz
cd - || exit

# Calculate the abundance of all organisms in the fastq file
sourmash gather \
    ${SOURMASH_PARAMS}\
    -k ${hash_kmer} \
    -o ${SOURMASH_OUTPUT}/${METHOD}/${SAMPLE_NAME}.profile.csv \
    ${SOURMASH_OUTPUT}/compute/${SAMPLE_NAME}-reads.sig \
    ${LOCAL_DB_PATH}/${LOCAL_DBNAME}.sbt.json |\
    tee -a ${LOG_DIR}/sourmash_${METHOD}.log

# Can also try for each genome only
# https://sourmash.readthedocs.io/en/latest/tutorial-basic.html?#compare-reads-to-assemblies

aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
