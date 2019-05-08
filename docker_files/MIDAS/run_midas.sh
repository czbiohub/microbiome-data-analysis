#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

####################################################################
#
# Requirements
#   vcpus: 4
#   memory: 20GB
#   volume:
#        default DB: 36GB
#
# Expected Variables:
#   Required:
#       coreNum=4
#       fastq1="s3://..."
#       s3OutputPath="s3://..."
#
#   Optional (can be left blank):
#       fastq2="s3://..." OR ""
#       hardTrim=80 OR ""
#       subsetReads=4000000 OR ""
#       s3path2db="s3://..." OR empty -> DEFAULT: "s3://czbiohub-brianyu/Synthetic_Community/Genome_References/MIDAS/1.2/midas_db_v1.2.tar.gz"
#
#   Sample Batch Submission:
#   aegea batch submit --queue aegea_batch --image sunitjain/midas:latest \
#     --storage /mnt=500 --memory 64000 --vcpus 16 \
#     --command="export coreNum=16; \
#     export fastq1=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz; \
#     export fastq2=s3://czbiohub-brianyu/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz; \
#     export s3OutputPath=s3://czbiohub-brianyu/Sunit_Jain/20190304_Midas_Test/Dorea-longicatena-DSM-13814/; \
#     export subsetReads=''; \
#     export hardTrim=''; \
#     ./run_midas.sh"
####################################################################

set -e
set -u
set -o pipefail

############################# SETUP ################################
START_TIME=$SECONDS
# export $@
export PATH="/opt/conda/bin:${PATH}"
export PYTHONPATH="/opt/conda/lib/python2.7/site-packages/midas"

echo $PATH
LOCAL=$(pwd)

OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
REF_DB="${OUTPUTDIR}/reference"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
QC_FASTQ="${LOCAL_OUTPUT}/trimmed_fastq"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
defaultDB="s3://czbiohub-microbiome/ReferenceDBs/Midas/v1.2/midas_db_v1.2.tar.gz"

# remove trailing '/'
s3OutputPath=${s3OutputPath%/}
SAMPLE_NAME=$(basename ${s3OutputPath})
SPECIES_OUT="${LOCAL_OUTPUT}/midas/${SAMPLE_NAME}"
GENES_OUT="${LOCAL_OUTPUT}/midas/${SAMPLE_NAME}"

mkdir -p ${OUTPUTDIR} ${RAW_FASTQ} ${QC_FASTQ} ${SPECIES_OUT} ${GENES_OUT} ${LOG_DIR}

trap '{ 
    aws s3 sync ${LOCAL_OUTPUT}/ ${s3OutputPath}/;
    rm -rf ${OUTPUTDIR} ; 
    exit 255; 
    }' 1 

############################ DATABASE ##############################

## Database
if [[ -z ${s3path2db+defaultDB} ]]; then
    echo "Missing s3 path to database. Using default."
    s3path2db=${defaultDB}
fi

localPath2db_tgz=${REF_DB}/$(basename ${s3path2db})
localPath2db=${REF_DB}/$(basename ${s3path2db} .tar.gz)
aws s3 cp --quiet ${s3path2db} ${localPath2db_tgz}

# Decompress
tar -xzf ${localPath2db_tgz} --directory ${REF_DB}/
rm -rf ${localPath2db_tgz}

########################## PARAMETERS ##############################

# Constant definitions for bbduk
trimQuality=${trimQuality:-25}
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

## Species
MIDAS_SPECIES_PARAMS="-t ${coreNum} -d ${localPath2db} --remove_temp"
MIDAS_GENES_PARAMS=${MIDAS_SPECIES_PARAMS}

# Use bbduk to trim reads, -eoom exits when out of memory
BBDUK_PARAMS="ref=adapters ktrim=r k=${kmer_value} mink=${min_kmer_value} hdist=1 tbo qtrim=rl trimq=${trimQuality} minlen=${minLength}"

### Single end
if [[ -n ${fastq1} ]]; then
    FASTQ1_FILENAME=$(basename ${fastq1} .gz)
    
    BBDUK_PARAMS="${BBDUK_PARAMS} in1=${RAW_FASTQ}/${FASTQ1_FILENAME} out1=${QC_FASTQ}/${FASTQ1_FILENAME}"
    MIDAS_SPECIES_PARAMS="${MIDAS_SPECIES_PARAMS} -1 ${QC_FASTQ}/${FASTQ1_FILENAME}"
    MIDAS_GENES_PARAMS=${MIDAS_SPECIES_PARAMS}
    # Download
    aws s3 cp ${fastq1} ${RAW_FASTQ}/
    gunzip ${RAW_FASTQ}/${FASTQ1_FILENAME}.gz
else
    echo "[FATAL] Missing fastq file."
fi

### Paired
if [[ -n ${fastq2:-} ]]; then
    FASTQ2_FILENAME=$(basename ${fastq2} .gz)

    BBDUK_PARAMS="${BBDUK_PARAMS} in2=${RAW_FASTQ}/${FASTQ2_FILENAME} out2=${QC_FASTQ}/${FASTQ2_FILENAME}"
    MIDAS_SPECIES_PARAMS="${MIDAS_SPECIES_PARAMS} -2 ${QC_FASTQ}/${FASTQ2_FILENAME}"
    MIDAS_GENES_PARAMS=${MIDAS_SPECIES_PARAMS}
    # Download
    aws s3 cp ${fastq2} ${RAW_FASTQ}/
    gunzip ${RAW_FASTQ}/${FASTQ2_FILENAME}.gz
fi

### Hard Trimming
if [[ -n ${hardTrim:-} ]]; then
    MIDAS_SPECIES_PARAMS="${MIDAS_SPECIES_PARAMS} --read_length ${hardTrim}"
fi

### Read subsets
if [[ -n ${subsetReads:-} ]]; then
    MIDAS_SPECIES_PARAMS="${MIDAS_SPECIES_PARAMS} -n ${subsetReads}"
fi

########################### EXECUTION ##############################

BBDUK="bbduk.sh -Xmx16g -eoom ${BBDUK_PARAMS} | tee -a ${LOG_DIR}/bbduk.log.txt"
MIDAS_SPECIES="run_midas.py species ${SPECIES_OUT} ${MIDAS_SPECIES_PARAMS} | tee -a ${LOG_DIR}/midas_species.log.txt"
MIDAS_GENES="run_midas.py genes ${GENES_OUT} ${MIDAS_GENES_PARAMS} | tee -a ${LOG_DIR}/midas_genes.log.txt"

############################# BBDUK ################################

if eval "${BBDUK}"; then
    echo "[$(date)] BBDUK complete."
    aws s3 sync ${LOCAL_OUTPUT}/ ${s3OutputPath}/
else
    echo "[$(date)] BBDUK failed."
    exit 1;
fi

####################### MIDAS - Species ############################

if eval "${MIDAS_SPECIES}"; then
    echo "[$(date)] MIDAS Species complete. Syncing output to ${s3OutputPath}"
    aws s3 sync ${LOCAL_OUTPUT}/ ${s3OutputPath}/
else
    echo "[$(date)] MIDAS Species failed."
    exit 1;
fi

######################## MIDAS - Genes #############################

if eval "${MIDAS_GENES}"; then
    echo "[$(date)] MIDAS Genes complete. Syncing output to ${s3OutputPath}"
    aws s3 sync ${LOCAL_OUTPUT}/ ${s3OutputPath}/
else
    echo "[$(date)] MIDAS Genes failed."
    exit 1;
fi

######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${OUTPUTDIR}/job.complete
echo "Live long and prosper" >> ${OUTPUTDIR}/job.complete
############################ PEACE! ################################
aws s3 cp ${OUTPUTDIR}/job.complete ${s3OutputPath}/