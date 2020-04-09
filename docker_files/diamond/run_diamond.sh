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
DUP_REMOVAL=${DUP_REMOVAL:-false}
SKIP_QC=${SKIP_QC:-true}
LOCAL_DB_PATH=${LOCAL}/databases
PROGRAM=${PROGRAM:-blastx}

# fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/Bowtie2_Test/Dorea-longicatena-DSM-13814
# S3DBPREFIX=s3://czbiohub-microbiome/Synthetic_Community/Genome_References/ncbi_fasta/Dorea-longicatena-DSM-13814-GCF_000154065.1_ASM15406v1

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
TMP_DMND_OUTPUT="${OUTPUTDIR}/diamond"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
DMND_OUTPUT="${LOCAL_OUTPUT}/diamond"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})
FASTQ_FILES=${OUTPUTDIR}/combined.fastq.gz

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}"
mkdir -p "${LOCAL_DB_PATH}" "${DMND_OUTPUT}" "${TMP_DMND_OUTPUT}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1 

hash_kmer=${hash_kmer:-51}

if [ "${SKIP_QC}" = false ];then
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
        refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt" &> \
        ${LOG_DIR}/bbduk.log
else
    # Copy QC-ed fastq.gz files from S3, only 2 files per sample
    aws s3 cp --quiet ${fastq1} "${QC_FASTQ}/read1_trimmed.fastq.gz"
    aws s3 cp --quiet ${fastq2} "${QC_FASTQ}/read2_trimmed.fastq.gz"
fi

cat \
    ${QC_FASTQ}/read1_trimmed.fastq.gz \
    ${QC_FASTQ}/read2_trimmed.fastq.gz > \
    ${FASTQ_FILES}

##
LOCAL_DB_NAME=$(basename ${S3DBPATH})
# S3DBPATHDIR=$(dirname ${S3DBPATH})
aws s3 cp ${S3DBPATH} ${LOCAL_DB_PATH}/${LOCAL_DB_NAME}

OUTPUT_PREFIX="${SAMPLE_NAME}_vs_db_${LOCAL_DB_NAME}"
DIAMOND_LOG="${OUTPUT_PREFIX}.log"

## Map the reads
diamond \
    ${PROGRAM} \
    --db ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} \
    --dbsize 1000000 \
    --tmpdir ${TMP_DMND_OUTPUT} \
    --threads ${CORE_NUM} \
    --outfmt 6 \
    --out ${DMND_OUTPUT}/${OUTPUT_PREFIX}.${PROGRAM}.txt \
    --query ${FASTQ_FILES} \
    --log &> \
    ${LOG_DIR}/${DIAMOND_LOG}

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
rm -rf "${OUTPUTDIR}"