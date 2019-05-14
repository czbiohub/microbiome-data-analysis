#!/bin/bash -x

# aws s3 ls s3://czbiohub-microbiome/Synthetic_Community/2019.04.17_strains_only_accuAlign_analysis/ | awk -v s3path="s3://czbiohub-microbiome/Synthetic_Community/2019.04.17_strains_only_accuAlign_analysis/" '{printf "%s%s\n",s3path,$2}' > outputs_folder.list
S3PATH=${1}
S3PATH=${S3PATH%/}
SAMPLE_FOLDER_NAME=$(basename ${S3PATH})
aws s3 cp ${S3PATH}/accu_align/ ${SAMPLE_FOLDER_NAME}/ --recursive --exclude '*' --include '*processed*sortedByCoord.bam'
for bamFile in $(find ${SAMPLE_FOLDER_NAME}/ -name "*bam"); do
    BAM_BASE=$(basename ${bamFile} .bam)
    bedtools genomecov -ibam ${bamFile} 1> ${SAMPLE_FOLDER_NAME}/${BAM_BASE}.genomecov.bed.txt
    aws s3 cp ${SAMPLE_FOLDER_NAME}/${BAM_BASE}.genomecov.bed.txt ${S3PATH}/genomeCoverage/
done
rm -rf ${SAMPLE_FOLDER_NAME}/