#!/bin/bash -x

set -e
set -u

SAMPLE_NAME="${1}"
INPUT_S3PATH="s3://czbiohub-microbiome/Synthetic_Community/2019.04.17_strains_only_accuAlign_analysis"
OUTPUT_S3PATH="s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/bowtie_accu_align/2019-04-30_Alignment_tuning"
DIST_SCRIPT="./get_aln_distribution.py"

NUM_CORES=8

mkdir -p ${SAMPLE_NAME}

aws s3 cp --quiet ${INPUT_S3PATH}/${SAMPLE_NAME}/accu_align/${SAMPLE_NAME}.processed.bam ./${SAMPLE_NAME}/

samtools sort -@ ${NUM_CORES} -o ${SAMPLE_NAME}/${SAMPLE_NAME}.processed.sortedByCoord.bam ${SAMPLE_NAME}/${SAMPLE_NAME}.processed.bam
samtools index ${SAMPLE_NAME}/${SAMPLE_NAME}.processed.sortedByCoord.bam

python ${DIST_SCRIPT} ${SAMPLE_NAME}/${SAMPLE_NAME}.processed.sortedByCoord.bam ${SAMPLE_NAME}/${SAMPLE_NAME}.processed.sortedByCoord.csv

rm -rf ${SAMPLE_NAME}/${SAMPLE_NAME}.processed.bam ${SAMPLE_NAME}/${SAMPLE_NAME}.processed.sortedByCoord.bam ${SAMPLE_NAME}/${SAMPLE_NAME}.processed.sortedByCoord.bam.bai

aws s3 sync --quiet ./${SAMPLE_NAME}/ ${OUTPUT_S3PATH}/${SAMPLE_NAME}/