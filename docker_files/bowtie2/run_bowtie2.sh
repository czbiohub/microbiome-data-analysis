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
S3DBFASTA=${S3DBFASTA:-''}
S3DBPREFIX=${S3DBPREFIX:-''}
# Bowtie2 Parameters
maxInsert=${maxInsert:-3000}
maxAlignments=${maxAlignments:-300}

# S3DBFASTA=s3://czbiohub-microbiome/Synthetic_Community/Genome_References/ncbi_fasta/Dorea-longicatena-DSM-13814-GCF_000154065.1_ASM15406v1.fna
# fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/Bowtie2_Test/Dorea-longicatena-DSM-13814
# S3DBPREFIX=s3://czbiohub-microbiome/Synthetic_Community/Genome_References/ncbi_fasta/Dorea-longicatena-DSM-13814-GCF_000154065.1_ASM15406v1

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
TMP_BWT_OUTPUT="${OUTPUTDIR}/bowtie2"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
BWT_OUTPUT="${LOCAL_OUTPUT}/bowtie2"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}"
mkdir -p "${LOCAL_DB_PATH}" "${BWT_OUTPUT}" "${TMP_BWT_OUTPUT}"
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

##
if [ -n "${S3DBFASTA}" ]; then
    LOCAL_REFSEQ_NAME=$(basename ${S3DBFASTA})
    LOCAL_REFSEQ_EXT=${S3DBFASTA##*.}
    LOCAL_DB_NAME=$(basename ${S3DBFASTA} .${LOCAL_REFSEQ_EXT})
    aws s3 cp ${S3DBFASTA} ${LOCAL_DB_PATH}/

    ## Build the database
    bowtie2-build --seed 1712 \
        --threads ${coreNum} \
        ${LOCAL_DB_PATH}/${LOCAL_REFSEQ_NAME} \
        ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} |\
        tee -a ${LOG_DIR}/bowtie2_build_db_index.log
elif [ -n "${S3DBPREFIX}" ]; then
    LOCAL_DB_NAME=$(basename ${S3DBPREFIX})
    S3DBPATH=$(dirname ${S3DBPREFIX})
    aws s3 sync ${S3DBPATH} ${LOCAL_DB_PATH}
fi

OUTPUT_PREFIX="${SAMPLE_NAME}_vs_db_${LOCAL_DB_NAME}"

## Map the reads
# -t -D 10 -R 2 -L 31 -i S,0,2.50 -N 0 \
bowtie2 \
    --seed 1712 \
    --align-paired-reads \
    --very-sensitive \
    -X ${maxInsert} \
    -k ${maxAlignments} \
    --threads ${coreNum} \
    -x ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} \
    --no-mixed \
    --no-discordant \
    --end-to-end \
    --no-unal \
    --un-gz "${QC_FASTQ}/unpaired_unaligned.fastq.gz" \
    --un-conc-gz "${QC_FASTQ}/paired_discordant.fastq.gz" \
    -1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
    -2 "${QC_FASTQ}/read2_trimmed.fastq.gz" | \
samtools view \
    -bh \
    -@ ${coreNum} \
    -o ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.bam - | \
    tee -a ${LOG_DIR}/read_mapping.log

# Remove PCR duplicates
if [ "${DUP_REMOVAL}" = true ]; then
    samtools sort -n -o ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.bam -O BAM -@ ${coreNum} ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.bam
    samtools fixmate -cm ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.bam ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.fixmate.bam
    samtools sort -o ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.fixmate.bam -O BAM -@ ${coreNum} ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.fixmate.bam
    samtools markdup -s -S ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.fixmate.bam ${BWT_OUTPUT}/${OUTPUT_PREFIX}.name_sorted.markdup.bam
fi

samtools sort -o ${BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.bam -O BAM -@ ${coreNum} ${TMP_BWT_OUTPUT}/${OUTPUT_PREFIX}.bam
samtools index -@ ${coreNum} ${BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.bam ${BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.bam.bai

# Calculate alignment stats
# The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads.
samtools idxstats -@ ${coreNum} ${BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.bam &> ${BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.idxstats.txt
samtools stats -@ ${coreNum} ${BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.bam &> ${BWT_OUTPUT}/${OUTPUT_PREFIX}.coord_sorted.stats.txt

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
rm -rf "${OUTPUTDIR}"
