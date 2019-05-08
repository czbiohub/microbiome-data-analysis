#!/bin/bash -x
# USAGE: bash -x create_seedfile_for_submission.sh s3://czbiohub-microbiome/Original_Sequencing_Data/181203_A00111_0236_BHCWTNDSXX/Kazuki_Nagashima/

S3PATH="${1%/}"
HEADER="sampleName,fastq1,fastq2"

aws s3 ls "${S3PATH}/" | sort -u | awk '{print $4}' | sort > $$.check.txt
aws s3 ls "${S3PATH}/" | sort -u | awk -v s3path="${S3PATH}" '{printf "%s/%s\n",s3path,$4}' | sort > $$.s3paths.txt

echo "${HEADER}" > seedfile.txt
paste $$.check.txt $$.s3paths.txt | paste - - | cut -f 1,2,4 | sed -E "s/_S[0-9]+_R1_001.fastq.gz//" | awk '{printf "%s,%s,%s\n",$1,$2,$3}' >> seedfile.txt

rm $$.check.txt $$.s3paths.txt