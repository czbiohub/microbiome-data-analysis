#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-15};
LOCAL_DB_PATH=${LOCAL}/databases
METHOD=${METHOD:-"search"}

S3DBPATH=s3://czbiohub-microbiome/Synthetic_Community/Genome_References/ncbi_fasta/Dorea-longicatena-DSM-13814-GCF_000154065.1_ASM15406v1.fna
fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/Sourmash_Test/Dorea-longicatena-DSM-13814

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
SOURMASH_OUTPUT="${LOCAL_OUTPUT}/sourmash"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

mkdir -p "${OUTPUTDIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${LOCAL_OUTPUT}" "${SOURMASH_OUTPUT}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1 

# Get database
mkdir -p "${LOCAL_DB_PATH}"
cd "${LOCAL_DB_PATH}" || exit

if [[ ${METHOD} == "search" ]]; then
    if [[ -n ${S3DBPATH:-} ]]; then
        LOCAL_DBNAME=$(basename ${S3DBPATH} .gz)
        mkdir -p ${SOURMASH_OUTPUT}/compute
        mkdir -p ${SOURMASH_OUTPUT}/search
        LOCAL_SM_DB="${SOURMASH_OUTPUT}/compute/${LOCAL_DBNAME}-ref.sig"
        aws s3 --quiet cp ${S3DBPATH} .
        # gunzip "${LOCAL_DBNAME}.gz"
        # rm -rf "${LOCAL_DBNAME}.gz"
        # Compute the signature of the input fastq file
        sourmash compute --scaled 1000 -k 51 ${LOCAL_DBNAME} \
            --merge ${LOCAL_DBNAME} -o ${LOCAL_SM_DB}

        SOURMASH_PARAMS="--containment"
    else
        echo "[FATAL] Using 'sourmash search', but missing reference genome to search."
        exit 1
    fi
elif [[ ${METHOD} = "gather" ]]; then
    # DEFAULT_REF_S3DB="s3://czbiohub-microbiome/ReferenceDBs/Sourmash/Genbank/2018-03/genbank-d2-k51.tar.gz"
    # S3DBPATH=${S3DBPATH:-$DEFAULT_REF_S3DB}
    # LOCAL_DBNAME=$(basename ${S3DBPATH} .tar.gz)
    # LOCAL_SM_DB="${SOURMASH_OUTPUT}/compute/${LOCAL_DBNAME}-ref.sig"
    # aws s3 cp ${S3DBPATH} .
    
    # TODO: Needs EFS.
    echo "[FATAL] This section of the code is under construction. Please use the default method 'search' for now."
    exit 1;
fi

cd "${LOCAL}" || exit

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

# Use khmer to remove low abundance reads that may be artifacts (optional step)
trim-low-abund.py -C 3 -Z 18 -V -M 10e9 \
    --summary-info tsv \
    "${QC_FASTQ}/read1_trimmed.fastq.gz" \
    "${QC_FASTQ}/read2_trimmed.fastq.gz"

# Compute the signature of the input fastq file
sourmash compute --scaled 1000 -k 51 \
    read*_trimmed*.abundtrim \
    --merge ${SAMPLE_NAME} \
    -o ${SOURMASH_OUTPUT}/compute/${SAMPLE_NAME}-reads.sig

rm -rf read*_trimmed*.abundtrim

# Calculate the abundance of all organisms in the fastq file
sourmash ${METHOD} \
    ${SOURMASH_PARAMS}\
    -k 51 \
    -o ${SOURMASH_OUTPUT}/${METHOD}/${SAMPLE_NAME}.profile.csv \
    ${LOCAL_SM_DB} \
    ${SOURMASH_OUTPUT}/compute/${SAMPLE_NAME}-reads.sig

# Can also try for each genome only
# https://sourmash.readthedocs.io/en/latest/tutorial-basic.html?#compare-reads-to-assemblies

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
aws s3 sync --quiet "${LOCAL_OUTPUT}" ${S3OUTPUTPATH}
# rm -rf "${OUTPUTDIR}"
