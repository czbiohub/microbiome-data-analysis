#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u
set -o pipefail

####################################################################
# accu_align.sh
#
# Accurate alignment of reads from synthetic microbiome to a comebined
# set of microbial genomes using bowtie2. The accuracy comes from taking
# care of reads that map equally well to multiple references.
#
# This script is used in conjunction with aegea batch
#
# Variables required from sourcing script
# coreNum=4; numPerCore=1G; maxInsert=3000; maxAlignments=200;
# S3DBPATH=/czbiohub-microbiome/Synthetic_Community/Genome_References/Bowtie2Index_090718
# SAMPLE_NAME; fastq1=/czbiohub-microbiome/...; fastq2;
# BAM_OUTPUT; REL_AB_OUTPUT; READ_ACC_OUTPUT
#
# bbtools assumes 16G of memory -Xmx16g; needs sambamba in conda env
#
# Revision History:
# 2018.09.08 Brian Yu Created
# 2018.09.09 Brian Yu Added alignment organization and tabulate_alignment_fragment.py
# 2018.12.09 Brian Yu Added off limit region removal, coverage and fragment tabulate
#####################################################################

START_TIME=$SECONDS
# export $@
export PATH="/opt/conda/bin:${PATH}"

coreNum="${coreNum:-15}"
memPerCore="${memPerCore:-2G}"
maxInsert="${maxInsert:-3000}"
maxAlignments="${maxAlignments:-200}"
minPercId="${minPercId:-0}"
minReadQuality="${minReadQuality:-0}"
minMapQuality="${minMapQuality:-10}"
minAlnCov="${minAlnCov:-0}"

# Inputs
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/bowtie_accu_align/2019-04-30_StrainVerification/Dorea-longicatena-DSM-13814
# fastq1=s3://czbiohub-microbiome/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-microbiome/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz

S3DBPATH=s3://czbiohub-microbiome/Synthetic_Community/Genome_References/Bowtie2Index_090718
REFDBNAME=combined_104_reference_genomes

SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

echo $PATH
LOCAL=$(pwd)

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
TMP_OUTPUTS="${OUTPUTDIR}/bowtie2"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
STATS_DIR="${LOCAL_OUTPUT}/Stats"
ACCU_ALN_OUTPUT="${LOCAL_OUTPUT}/accu_align"
GENOME_COV_OUTPUT="${LOCAL_OUTPUT}/genome_coverage"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
S3DBPATH=${S3DBPATH%/}
LOCAL_DB_PATH="${OUTPUTDIR}/reference"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${STATS_DIR}"
mkdir -p "${LOCAL_DB_PATH}" "${ACCU_ALN_OUTPUT}" "${TMP_OUTPUTS}" "${GENOME_COV_OUTPUT}"

trap '{ aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}";
    rm -rf ${OUTPUTDIR} ; 
    exit 255; }' 1 

adapterFile="adapters,phix"
# offLimitRegions="./data/combined_excluded_regions_threshold9.bed"
scriptFolder="./scripts"
DBNAME=${LOCAL_DB_PATH}/${REFDBNAME}

# Copy genome reference over
aws s3 sync --quiet ${S3DBPATH}/ ${LOCAL_DB_PATH}/
referenceNameFile=${LOCAL_DB_PATH}/genome_name_list.csv

# Constant definitions for bbduk
trimQuality="${trimQuality:-25}"
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

echo "Starting to Process Sample: "${SAMPLE_NAME}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} ${RAW_FASTQ}/read1.fastq.gz
aws s3 cp --quiet ${fastq2} ${RAW_FASTQ}/read2.fastq.gz

# Use bbduk to trim reads, -eoom exits when out of memory
bbduk.sh -Xmx16g -eoom \
    in1=${RAW_FASTQ}/read1.fastq.gz \
    in2=${RAW_FASTQ}/read2.fastq.gz \
    out1=${QC_FASTQ}/read1_trimmed.fastq.gz \
    out2=${QC_FASTQ}/read2_trimmed.fastq.gz \
    ref=${adapterFile} \
    ktrim=r \
    k=${kmer_value} \
    mink=${min_kmer_value} \
    hdist=1 tbo qtrim=rl \
    trimq=${trimQuality} \
    minlen=${minLength} \
    refstats=${STATS_DIR}/adapter_trimming_stats_per_ref.txt |\
    tee -a ${LOG_DIR}/bbduk.log.txt
	
# bowtie2 alignment returning multiple alignments and using longer max insert size limites
# output samtools bam file with only properly aligned paired reads.
bowtie2 -t -D 10 -R 2 -L 31 -i S,0,2.50 -N 0 \
    -X ${maxInsert} \
    -k ${maxAlignments} \
    --threads ${coreNum} \
    -x ${DBNAME} \
    --no-mixed \
    --no-discordant \
    --end-to-end \
    --no-unal \
    -1 ${QC_FASTQ}/read1_trimmed.fastq.gz \
    -2 ${QC_FASTQ}/read2_trimmed.fastq.gz | \
    samtools view \
        -@ ${coreNum} \
        -bh \
        -f 3 \
        -o ${ACCU_ALN_OUTPUT}/${SAMPLE_NAME}.bam - |\
    tee -a ${LOG_DIR}/read_mapping.log.txt

# Mark PCR Duplicates
samtools sort \
    -n \
    -@ ${coreNum} \
    -m ${memPerCore} \
    -O BAM \
    ${ACCU_ALN_OUTPUT}/${SAMPLE_NAME}.bam |\
samtools fixmate \
    -O BAM \
    -cm \
    - - | \
samtools sort \
    -@ ${coreNum} \
    -m ${memPerCore} \
    -O BAM \
    - |\
samtools markdup \
    -s \
    -S - \
    ${TMP_OUTPUTS}/${SAMPLE_NAME}.coord_sorted.markdup.bam |\
    tee -a ${LOG_DIR}/samtools_markdup.log.txt

# 3328 =
#   not primary alignment (0x100)
#   read is PCR or optical duplicate (0x400)
#   supplementary alignment (0x800)
# samtools view -F 3328 -q 10 Dorea-longicatena-DSM-13814.processed.bam | cut -f1 | sort | uniq | wc -l

# Remove reads that fall into off limit regions
# bedtools intersect \
#     -abam ${TMP_OUTPUTS}/${SAMPLE_NAME}.coord_sorted.markdup.bam \
#     -v -b ${offLimitRegions} 1> \
#     ${TMP_OUTPUTS}/${SAMPLE_NAME}.coord_sorted.markdup.offlim_rm.bam |\
#     tee -a ${LOG_DIR}/bedtools_offlimits.log.txt

OUTPUT_PREFIX="${SAMPLE_NAME}.processed"
    
samtools sort \
    -n \
    -@ ${coreNum} \
    -m ${memPerCore} \
    -O BAM \
    -o ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.bam \
    ${TMP_OUTPUTS}/${SAMPLE_NAME}.coord_sorted.markdup.bam |\
    tee -a ${LOG_DIR}/samtools_sort_by_name_2.log.txt
# ${TMP_OUTPUTS}/${SAMPLE_NAME}.coord_sorted.markdup.offlim_rm.bam |\

# Split bam files into 3 bam files
python ${scriptFolder}/split_multiple_alignment.py \
    -bam ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.bam \
    -max_align ${maxAlignments} \
    -min_pid ${minPercId} \
    -min_readq ${minReadQuality} \
    -min_mapq ${minMapQuality} \
    -min_aln_cov ${minAlnCov} |\
    tee -a ${LOG_DIR}/split_multiple_alignment.log.txt

# Tabulate read count
totalReads=$(( $( zcat ${RAW_FASTQ}/read1.fastq.gz | wc -l ) / 4 ))
readsAfterTrim=$(( $( zcat ${QC_FASTQ}/read1_trimmed.fastq.gz | wc -l ) / 4 ))
uniqueReads=$( samtools view -f 0x40 ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.unique_alignments.bam | cut -f1 | sort -u | wc -l )
multipleReads1=$( samtools view -f 0x40 ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.multiple_alignments_multiple_genomes.bam | cut -f1 | sort -u | wc -l )
multipleReads2=$( samtools view -f 0x40 ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.multiple_alignments_unique_genome.bam | cut -f1 | sort -u | wc -l )
multipleReads=$(( multipleReads1 + multipleReads2 ))
echo 'Sample_Name,Total_Fragments,Fragments_After_Trim,Fragments_Aligned_Uniquely,Fragments_Aligned_Multiple_Times' > ${ACCU_ALN_OUTPUT}/read_accounting.csv
echo ${SAMPLE_NAME}','${totalReads}','${readsAfterTrim}','${uniqueReads}','${multipleReads} >> ${ACCU_ALN_OUTPUT}/read_accounting.csv

# Go through the multimapped reads to assign them using a python script.
date
# bedtools genomecov \
#     -ibam ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.unique_alignments.bam 1> \
#     ${GENOME_COV_OUTPUT}/${OUTPUT_PREFIX}.unique_alignments.bed |\
#     tee -a ${LOG_DIR}/genomeCov_unique_alignment.log.txt

# # multiple 1
# bedtools genomecov \
#     -ibam ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.multiple_alignments_multiple_genomes.bam 1> \
#     ${GENOME_COV_OUTPUT}/${OUTPUT_PREFIX}.multiple_alignments_multiple_genomes.bed |\
#     tee -a ${LOG_DIR}/genomeCov_multiple_alignments_multiple_genomes.log.txt

# # multiple 2
# bedtools genomecov \
#     -ibam ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.multiple_alignments_unique_genome.bam 1> \
#     ${GENOME_COV_OUTPUT}/${OUTPUT_PREFIX}.multiple_alignments_unique_genome.bed |\
#     tee -a ${LOG_DIR}/genomeCov_multiple_alignments_unique_genome.log.txt

python ${scriptFolder}/tabulate_alignment_fragment.py \
    -s ${SAMPLE_NAME} \
    ${DBNAME}.fasta \
    ${referenceNameFile} \
    ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.unique_alignments.bam \
    ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.multiple_alignments_multiple_genomes.bam \
    ${ACCU_ALN_OUTPUT}/${OUTPUT_PREFIX}.multiple_alignments_unique_genome.bam \
    ${ACCU_ALN_OUTPUT}/tabulated_alignment_fragment.csv
date

pwd
echo "Alignment completed."
ls ${LOCAL}
du -sh ${LOCAL}
date
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"