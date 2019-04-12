#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u

START_TIME=$SECONDS

LOCAL=$(pwd)
coreNum=${coreNum:-15};
IGG_CODE_PATH=${LOCAL}/src
IGG_DB_PATH=${LOCAL}/databases
DEFAULT_REF_DB="s3://czbiohub-microbiome/ReferenceDBs/IGGdb/gut_only/v1.0.0/iggdb_v1.0.0_gut.tar.gz"
S3DBPATH=${S3DBPATH:-$DEFAULT_REF_DB}
LOCAL_DBNAME=$(basename ${S3DBPATH} .tar.gz)
# sampleName=Dorea-longicatena-DSM-13814
# fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
# s3OutputPath=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/IGGsearch_Test/Dorea-longicatena-DSM-13814

# Get Code
wget -q "https://github.com/snayfach/IGGsearch/archive/v1.0.0.tar.gz"
mkdir -p "${IGG_CODE_PATH}"
tar xzf v1.0.0.tar.gz --directory "${IGG_CODE_PATH}"

export PYTHONPATH="${IGG_CODE_PATH}/IGGsearch-1.0.0/iggsearch"
export PATH="${IGG_CODE_PATH}/IGGsearch-1.0.0:/opt/conda/bin:${PATH}"

# Get database
mkdir -p "${IGG_DB_PATH}"
cd "${IGG_DB_PATH}" || exit
aws s3 cp --quiet ${S3DBPATH} .
tar xzf "${LOCAL_DBNAME}.tar.gz" --directory "${IGG_DB_PATH}"
cd "${LOCAL}" || exit

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
s3OutputPath=${s3OutputPath%/}
SAMPLE_NAME=$(basename ${s3OutputPath})

mkdir -p "${OUTPUTDIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${LOCAL_OUTPUT}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1 

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
    refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt"

run_iggsearch.py search \
    --m1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
    --m2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
    --db_dir "${IGG_DB_PATH}/${LOCAL_DBNAME}" \
    --outdir "${LOCAL_OUTPUT}/${SAMPLE_NAME}" \
    --threads "${coreNum}"

# # To roll up counts to other taxonomic ranks
# run_iggsearch.py reformat \
#     --indir "${LOCAL_OUTPUT}/IGG_search" \
#     --outdir "${LOCAL_OUTPUT}/IGG_reformat" \
#     --db_dir "${IGG_DB_PATH}/${LOCAL_DBNAME}" \
#     --taxdb gtdb \
#     --taxrank family
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
aws s3 sync --quiet "${LOCAL_OUTPUT}" ${s3OutputPath}
# rm -rf "${OUTPUTDIR}"
