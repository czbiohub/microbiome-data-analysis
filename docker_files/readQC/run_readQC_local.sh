#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u
set -o pipefail

export PATH="/opt/conda/bin:${PATH}"

START_TIME=$SECONDS


LOCAL=$(pwd)
coreNum=${coreNum:-15};

sampleName=${1}
fastq1=${2}
fastq2=${3:-""}

# sampleName=Dorea-longicatena-DSM-13814
# fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/IGGsearch_Test/Dorea-longicatena-DSM-13814
# Setup directory structure
# sampleName=${sampleName}
OUTPUTDIR=${LOCAL}/${sampleName}
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"

PRE_FASTQC="${LOCAL_OUTPUT}/preQC_stats"
QC_FASTQ="${LOCAL_OUTPUT}/trimmed_fastq"
POST_FASTQC="${LOCAL_OUTPUT}/postQC_stats"
# S3OUTPUTPATH=${S3OUTPUTPATH%/}
# SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

mkdir -p "${OUTPUTDIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${PRE_FASTQC}" "${POST_FASTQC}" "${LOCAL_OUTPUT}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1 

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
FASTQC_FILES="${RAW_FASTQ}/read1.fastq.gz"
TRIMMED_FASTQ="${QC_FASTQ}/read1_trimmed.fastq.gz"
BBDUK_INPUT="in1=${RAW_FASTQ}/read1.fastq.gz"
BBDUK_OUTPUT="out1=${QC_FASTQ}/read1_trimmed.fastq.gz"

if [[ -n ${fastq2} ]]; then
        aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"
        FASTQC_FILES="${FASTQC_FILES} ${RAW_FASTQ}/read2.fastq.gz"
        TRIMMED_FASTQ="${TRIMMED_FASTQ} ${QC_FASTQ}/read2_trimmed.fastq.gz"
        BBDUK_INPUT="${BBDUK_INPUT} in2=${RAW_FASTQ}/read2.fastq.gz"
        BBDUK_OUTPUT="${BBDUK_OUTPUT} out2=${QC_FASTQ}/read2_trimmed.fastq.gz"
fi

# PRE FASTQC
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd)":"$(pwd)" \
    quay.io/biocontainers/fastqc:0.11.9--0 \
        fastqc -o ${PRE_FASTQC} \
                -t ${coreNum} \
                --extract \
                ${FASTQC_FILES}

# Constant definitions for bbduk
adapterFile="adapters,phix"
trimQuality=${trimQuality:-25}
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

# Use bbduk to trim reads, -eoom exits when out of memory
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd)":"$(pwd)" \
        quay.io/biocontainers/bbmap:38.86--h1296035_0 \
        bbduk.sh -Xmx8g tbo -eoom hdist=1 qtrim=rl ktrim=r \
        ${BBDUK_INPUT} \
        ${BBDUK_OUTPUT} \
        ref=${adapterFile} \
        k="${kmer_value}" \
        mink="${min_kmer_value}" \
        trimq="${trimQuality}" \
        minlen="${minLength}" \
        refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt"

# POST FASTQC
docker container run --rm \
    --workdir "$(pwd)" \
    --volume "$(pwd)":"$(pwd)" \
    quay.io/biocontainers/fastqc:0.11.9--0 \
fastqc -o ${POST_FASTQC} \
        -t ${coreNum} \
        --extract \
        ${TRIMMED_FASTQ}

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
# aws s3 sync --quiet "${LOCAL_OUTPUT}" ${S3OUTPUTPATH}
# rm -rf "${OUTPUTDIR}"


