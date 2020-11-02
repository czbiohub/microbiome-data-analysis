#!/bin/bash -x

set -euo pipefail

THREADS=2
SAMPLE_NAME="${1}"
shift
GFF_FILE="${1}"
shift
BAM_FILES=( "$@" )

OUTPUT_DIR=${SAMPLE_NAME}/multicov
mkdir -p ${OUTPUT_DIR}

OUTPUT=${OUTPUT_DIR}/${SAMPLE_NAME}.multicov.tsv
GENES_GFF_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.genes.gff"

STRING_OF_BAMS=""
# foreach bam file:
for BAM_FILE in "${BAM_FILES[@]}"; do
    BAM_BASE_NAME=$(basename ${BAM_FILE} .coord_sorted.bam)

    FILTERED_BAM="${OUTPUT_DIR}/${BAM_BASE_NAME}.filtered.bam"
    SORTED_FILTERED_BAM="${OUTPUT_DIR}/${BAM_BASE_NAME}.filtered_coord_sorted.bam"

    ## filter all bams to 100% ID and aln_cov
    python filterBam.py -bam ${BAM_FILE} -out ${FILTERED_BAM}
    ## sort and index filtered bams
    samtools sort -@ ${THREADS} -T ${SAMPLE_NAME} -o ${SORTED_FILTERED_BAM} ${FILTERED_BAM}
    rm -rf ${FILTERED_BAM}
    samtools index -@ ${THREADS} ${SORTED_FILTERED_BAM}

    STRING_OF_BAMS="${STRING_OF_BAMS} ${SORTED_FILTERED_BAM}"
done

# filter gff file to only use the rows with 'genes' in the third column
zcat ${GFF_FILE} \
    | awk '$3 == "gene" {print $0}' &> ${GENES_GFF_FILE}

# compare reads from each filtered bam on gene features.
docker container run --rm \
    --workdir $(pwd) \
    --volume $(pwd):$(pwd) \
    quay.io/biocontainers/bedtools:2.26.0gx--he513fc3_4 \
    bedtools \
        multicov \
            -bams ${STRING_OF_BAMS} \
            -bed ${GENES_GFF_FILE} &> ${OUTPUT}

echo "All Done!"
echo "Huzzah!"