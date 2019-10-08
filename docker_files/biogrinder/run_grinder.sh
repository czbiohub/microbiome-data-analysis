#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
# S3FASTA="${1}"
# S3OUTPUTPATH="${2}"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
PIGZ_COMPRESSION_THREADS=${CORE_NUM:-4}

FASTA=$(basename -- "$S3FASTA")
PREFIX="${FASTA%.*}"

LOCAL=$(pwd)
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
RAW_FASTQ="${OUTPUTDIR}/interleaved_fastq"
PAIRED_FASTQ="${LOCAL_OUTPUT}/paired_fastq"
LOCAL_FASTA="${OUTPUTDIR}/${FASTA}"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${PAIRED_FASTQ}"
aws s3 cp ${S3FASTA} ${LOCAL_FASTA}

FWD="${PAIRED_FASTQ}/${PREFIX}.R1.fastq.gz"
REV="${PAIRED_FASTQ}/${PREFIX}.R2.fastq.gz"

# NUM_SEQS=$(grep -c ">" ${LOCAL_FASTA})
COV_FOLD=${COV_FOLD:-10}
MUT_DIST=${MUT_DIST:-"uniform 0"} # No errors
# USE: 
# for Illumina: MUT_DIST="poly4 3e-3 3.3e-8"
# for Sanger: MUT_DIST="linear 1 2"
# for Uniform: MUT_DIST="uniform 0.1"
# for No Errors: MUT_DIST="uniform 0" (?)

# Grinder Command
grinder \
    -cf ${COV_FOLD}\
    -rd 140 \
    -id 800 \
    -mo FR \
    -dc '-~*NX' \
    -md "${MUT_DIST}" \
    -rs 1712 \
    -am uniform \
    -ql 33 31 \
    -fq 1 \
    -od "${RAW_FASTQ}" \
    -bn "${PREFIX}" \
    -rf "${LOCAL_FASTA}" &> "${LOG_DIR}.log"

# Deinterleave and compress
paste - - - - - - - - < "${RAW_FASTQ}/${PREFIX}-reads.fastq" | sed "s/@/@${PREFIX}_/g" | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $FWD) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $REV

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > "${LOCAL_OUTPUT}/job.complete"
echo "Live long and prosper" >> "${LOCAL_OUTPUT}/job.complete"
############################ PEACE! ################################
## Sync output
aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"