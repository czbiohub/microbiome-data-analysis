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
####################################################################

set -e
set -u

############################# SETUP ################################

# export $@

LOCAL=$(pwd)

OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
SPECIES_OUT="${OUTPUTDIR}/midas_output"
GENES_OUT="${OUTPUTDIR}/midas_output"
REF_DB="${OUTPUTDIR}/reference"
# adapterFile="${LOCAL}/data/illumina_adapters.fa"
adapterFile="illumina_adapters.fa"
defaultDB="s3://czbiohub-brianyu/Synthetic_Community/Genome_References/MIDAS/1.2/midas_db_v1.2.tar.gz"
# remove trailing '/'

head ${adapterFile}

s3OutputPath=${s3OutputPath%/}
mkdir -p ${RAW_FASTQ} ${QC_FASTQ} ${SPECIES_OUT} ${GENES_OUT}

trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' EXIT 

######################### Downloads ################################

## Database
if [[ -z ${s3path2db+defaultDB} ]]; then
    echo "Missing s3 path to database. Using default."
    s3path2db=${defaultDB}
fi

localPath2db_tgz=${REF_DB}/$(basename ${s3path2db})
localPath2db=${REF_DB}/$(basename ${s3path2db} .tar.gz)

# Download
aws s3 cp --quiet ${s3path2db} ${localPath2db_tgz}

# Decompress
tar -xzf ${localPath2db_tgz} --directory ${REF_DB}/
rm -rf ${localPath2db_tgz}

# Download Fastq
FASTQ1_FILENAME=$(basename ${fastq1} .gz)
FASTQ2_FILENAME=$(basename ${fastq2} .gz)

aws s3 cp ${fastq1} ${RAW_FASTQ}/
gunzip ${RAW_FASTQ}/${FASTQ1_FILENAME}.gz

aws s3 cp ${fastq2} ${RAW_FASTQ}/
gunzip ${RAW_FASTQ}/${FASTQ2_FILENAME}.gz

########################## PARAMETERS ##############################

############################## QC ##################################

# Constant definitions for bbduk
trimQuality=25
minLength=50
kmer_value=23
min_kmer_value=11

####################### MIDAS - Species ############################

# Run MIDAS 
## Species
MIDAS_PARAMS="-t ${coreNum} -d ${localPath2db} --remove_temp"

# Use bbduk to trim reads, -eoom exits when out of memory
BBDUK_PARAMS="in1=${RAW_FASTQ}/${FASTQ1_FILENAME} out1=${QC_FASTQ}/${FASTQ1_FILENAME} ref=adapters ktrim=r k=${kmer_value} mink=${min_kmer_value} hdist=1 tbo qtrim=rl trimq=${trimQuality} minlen=${minLength}"

### Single end
if [[ -z ${fastq2} ]]; then
    MIDAS_PARAMS="${MIDAS_PARAMS} -1 ${QC_FASTQ}/${FASTQ1_FILENAME}"
else
### Paired
    MIDAS_PARAMS="${MIDAS_PARAMS} -1 ${QC_FASTQ}/${FASTQ1_FILENAME} -2 ${QC_FASTQ}/${FASTQ2_FILENAME}"
    BBDUK_PARAMS="${BBDUK_PARAMS} in2=${RAW_FASTQ}/${FASTQ2_FILENAME} out2=${QC_FASTQ}/${FASTQ2_FILENAME}"
fi

### Hard Trimming
if [[ -n ${hardTrim:-} ]]; then
    MIDAS_PARAMS="${MIDAS_PARAMS} --read_length ${hardTrim}"
fi

### Read subsets
if [[ -n ${subsetReads:-} ]]; then
    MIDAS_PARAMS="${MIDAS_PARAMS} -n ${subsetReads}"
fi

MIDAS="run_midas.py species ${SPECIES_OUT} ${MIDAS_PARAMS}"
BBDUK="bbduk.sh -Xmx16g -eoom ${BBDUK_PARAMS}"

########################### Execute ################################

if eval "${BBDUK}"; then
    echo "[$(date)] BBDUK complete."
else
    echo "[$(date)] BBDUK failed."
    exit 1;
fi

if eval "${MIDAS}"; then
    echo "[$(date)] MIDAS Species complete. Syncing output to ${s3OutputPath}/Species"
    aws s3 sync --quiet ${SPECIES_OUT}/species/ ${s3OutputPath}/Species/
else
    echo "[$(date)] MIDAS Species failed."
    exit 1;
fi
rm -rf ${SPECIES_OUT}
# TODO: ####################### MIDAS - Genes ##############################
