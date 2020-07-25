#!/bin/bash

####################################################################
# nanopore_basecall_demux.sh
#
# Variables required from sourcing script
# NANOPORE_RUNPATH, RUN_NAME
#
# Revision History:
# 2019.03.12 - created by Bryan Merrill
# 2019.08.28 - Dockerized
# 2020.02.04 - Adapting for basecall into czb-seqbot
# 2020.06.28 - Combined guppy_basecaller and guppy_barcoder calls. Streamlined code
# 2020.06.30 - Added --trim_barcodes --num_extra_bases_trim 2 for only barcoded demux
#              this is because guppy currently does not support adapter removal for non barcoded molecules
# 2020.07.02 - Add porechop functions and combine output fastq.gz
# 2020.07.03 - Added python script sample_sheet_process.py to do porechop.
#              Returns only 1 fastq.gz file for each barcode. If the barcode is in the sample sheet,
#              then the name of the file is Sample_Name.
#
#####################################################################

echo "**** START ****"
set -euo pipefail
START_TIME=$SECONDS
LOCAL=$(pwd)
NANOPORE_RUNPATH=${NANOPORE_RUNPATH:-s3://czb-seqbot/nanopore}
coreNum=${coreNum:-8}
guppy_threads=4
guppy_callers=$(( $coreNum / $guppy_threads ))

# Import samplesheet from S3
echo "Importing nanopore samplesheet ..."
# aws s3 cp s3://czb-seqbot/nanopore/nanopore-samplesheets/${RUN_NAME}.csv .
python sample_sheet_process.py --trial -t ${coreNum} -i ${RUN_NAME}.csv

# extract run parameters from sample sheet.
flowcell=$( grep -m1 "Flow Cell" ${RUN_NAME}.csv | cut -d, -f2 ) # ${flowcell:-FLO-MIN106}
kit=$( grep -m1 "Kit" ${RUN_NAME}.csv | cut -d, -f2 ) # ${kit:-SQK-LSK109}
barcoded=$( grep -m1 "Barcoded" ${RUN_NAME}.csv | cut -d, -f2 )
if [ $barcoded == "yes" ]; then
  barcodekit=$( grep -m1 "Barcode Kit" ${RUN_NAME}.csv | cut -d, -f2 ) # ${barcodekit:-EXP-NBD104 EXP-NBD114} # default is to find both kits
fi
guppy_threads=4
guppy_callers=$(( $coreNum / $guppy_threads ))

echo "Downloading nanopore fast5 files..."
aws s3 sync ${NANOPORE_RUNPATH}/${RUN_NAME}/ ${RUN_NAME}/ --quiet
echo "Download complete."

# Find all fast5 files and dump in new directory
mkdir -p raw_fast5
mv `find ${RUN_NAME}/ -name "*.fast5"` raw_fast5/

## Creating output directory
mkdir -p 00_MinIONQC
mkdir -p 01_BASECALLED
mkdir -p guppy_fastqs

## If not running with barcode
if [ $barcoded == "no" ]; then
  guppy_basecaller -i raw_fast5/ -s guppy_fastqs --flowcell $flowcell --kit $kit -x "cuda:0" --num_callers ${guppy_callers} --compress_fastq 2> 00_MinIONQC/barcoder.err 1> 00_MinIONQC/barcoder.out
  mkdir guppy_fastqs/barcode00
  mv guppy_fastqs/*.fastq.gz guppy_fastqs/barcode00/
fi
# If running with barcode
if [ $barcoded == "yes" ]; then
  guppy_basecaller -i raw_fast5/ -s guppy_fastqs --flowcell $flowcell --kit $kit -x "cuda:0" --num_callers ${guppy_callers} --trim_barcodes --num_extra_bases_trim 2 --barcode_kits "$barcodekit" --compress_fastq 2> 00_MinIONQC/barcoder.err 1> 00_MinIONQC/barcoder.out
fi

# Running QC on the output summary file
# -p $guppy_threads number of processors only useful if analyzing more than 1 summary file
Rscript MinIONQC.R -i guppy_fastqs/sequencing_summary.txt -o 00_MinIONQC -s FALSE
cp guppy_fastqs/sequencing_summary.txt 00_MinIONQC/
cp guppy_fastqs/sequencing_telemetry.js 01_BASECALLED/

# Change the files names to be the RUN_NAME
python sample_sheet_process.py -t ${coreNum} -i ${RUN_NAME}.csv

# Sync files back to S3
for folder in 00_MinIONQC 01_BASECALLED; do
  aws s3 sync --quiet $folder ${NANOPORE_RUNPATH}/${RUN_NAME}/$folder/
done

echo "File transfer complete."
pwd
echo "Nanopore basecalling and trim pipeline completed."
ls $LOCAL
du -sh $LOCAL
date
echo " **** END **** "
