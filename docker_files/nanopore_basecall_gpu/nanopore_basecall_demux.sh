#!/bin/bash

####################################################################
# nanopore_basecall_demux_wrapper.sh
#
#
# This script is used in conjunction with aegea batch (aegea_batch_basecall.sh).
#
# Variables required from sourcing script
# NANOPORE_RUNPATH, RUN_NAME, DATA_DEST, 
# coreNum, barcoded, flowcell, kit, barcodekit
#
# Revision History:
# 2019.03.12 - created by Bryan Merrill
# 2019.08.28 - Dockerized
# 
#####################################################################

echo "**** START ****"
set -euo pipefail
START_TIME=$SECONDS

LOCAL=$(pwd)
# The next two lines are necessary because there are two different /mnt locations
# Without this odd copy step, Snakemake fails (other things do too).
cp -pr * $LOCAL/
cd $LOCAL

export PATH="/opt/conda/bin:${PATH}"

coreNum=${coreNum:-16}
NANOPORE_RUNPATH=${NANOPORE_RUNPATH:-s3://czb-seqbot/nanopore} # UPDATE!!!
flowcell=${flowcell:-FLO-MIN106}
kit=${kit:-SQK-LSK109}
barcodekit=${barcodekit:-EXP-NBD104}
# Write something to remove trailing slashes from NANOPORE_RUNPATH, RUN_NAME, and DATA_DEST.


echo "Downloading nanopore fast5 files..."
aws s3 sync ${NANOPORE_RUNPATH}/${RUN_NAME}/ ${RUN_NAME}/ --quiet
echo "Download complete."

# Find all fast5 files and dump in new directory
mkdir -p raw_fast5
# Comment the line below if testing the script
mv `find ${RUN_NAME}/ -name "*.fast5"` raw_fast5/
# Uncomment line below if testing the script. Uncommenting the line will make things run fast enough to do debugging.
# mv `find ${RUN_NAME}/ -name "*.fast5" | head -n 2` raw_fast5/
# PROCESSING THE FILES
## Running guppy basecaller
guppy_threads=4
guppy_callers=$(( $coreNum / $guppy_threads ))
guppy_basecaller -i raw_fast5/ -s guppy_fastqs --flowcell $flowcell --kit $kit -x "cuda:0" #--cpu_threads_per_caller ${guppy_threads} --num_callers ${guppy_callers}

## Combining fastq files into one
mkdir fastq_concat
cat guppy_fastqs/*.fastq > fastq_concat/fastq_runid_000000000000_0.fastq

## Running QC on the output summary file
mkdir 00_MinIONQC
Rscript MinIONQC.R -i guppy_fastqs/sequencing_summary.txt -o 00_MinIONQC
cp -p guppy_fastqs/sequencing_summary.txt 00_MinIONQC/

## If run has barcodes, do demultiplexing
mkdir 01_BASECALLED
if [ $barcoded == "yes" ]; then
guppy_barcoder -r -i fastq_concat -s 01_BASECALLED -t $coreNum --barcode_kits $barcodekit 2> 00_MinIONQC/barcoder.err 1> 00_MinIONQC/barcoder.out
fi

if [ $barcoded == "no" ]; then
mv fastq_concat 01_BASECALLED/barcode00
fi

pigz 01_BASECALLED/barcode*/*.fastq

for dir in `ls -d 01_BASECALLED/barcode*/`; do 
barcode=`echo $dir | sed 's|01_BASECALLED/barcode||;s|/||'`
for file in `ls $dir*`; do 
newname=${RUN_NAME}__barcode${barcode}__`echo $file | awk -F'_' '{print $5}'`; mv $file $dir$newname
done
done


for folder in 00_MinIONQC 01_BASECALLED; do
aws s3 sync $folder ${DATA_DEST}/${RUN_NAME}/$folder/;
done

echo "File transfer complete."
pwd
echo "Nanopore basecalling pipeline completed."
ls $LOCAL
du -sh $LOCAL
date
echo " **** END **** "
