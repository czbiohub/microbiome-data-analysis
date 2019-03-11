#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u

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
# genomeReferenceDir=/czbiohub-brianyu/Synthetic_Community/Genome_References/Bowtie2Index_090718
# sampleName; fastq1=/czbiohub-brianyu/...; fastq2;
# bamOutput; relativeAbundanceOutput; readAccountingOutput
#
# bbtools assumes 16G of memory -Xmx16g; needs sambamba in conda env
#
# Revision History:
# 2018.09.08 Brian Yu Created
# 2018.09.09 Brian Yu Added alignment organization and tabulate_alignment_fragment.py
# 2018.12.09 Brian Yu Added off limit region removal, coverage and fragment tabulate
#####################################################################

# Prepare the EC2 instance with the ECR image
# /mnt has a volume of defined size mounted
LOCAL=/mnt

adapterFile=${LOCAL}/data/illumina_adapters.fa
offLimitRegions=${LOCAL}/data/combined_excluded_regions_threshold9.bed
scriptFolder=${LOCAL}/scripts

# Copy genome reference over
aws s3 sync --quiet s3:/${genomeReferenceDir}/ ${LOCAL}/reference/
bowtie2GenomeBase=${LOCAL}/reference/combined_104_reference_genomes
referenceNameFile=${LOCAL}/reference/genome_name_list.csv

# Constant definitions for bbduk
trimQuality=25
minLength=50
kmer_value=23
min_kmer_value=11

# Set {tempFolder} definition
tempFolder=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
mkdir ${tempFolder}

echo "Starting to Process Sample: "${sampleName}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet s3:/${fastq1} ${tempFolder}/read1.fastq.gz
aws s3 cp --quiet s3:/${fastq2} ${tempFolder}/read2.fastq.gz

# Use bbduk to trim reads, -eoom exits when out of memory
bbduk.sh -Xmx16g -eoom in1=${tempFolder}/read1.fastq.gz in2=${tempFolder}/read2.fastq.gz out1=${tempFolder}/read1_trimmed.fastq.gz out2=${tempFolder}/read2_trimmed.fastq.gz ref=${adapterFile} ktrim=r k=${kmer_value} mink=${min_kmer_value} hdist=1 tbo qtrim=rl trimq=${trimQuality} minlen=${minLength}
	
# bowtie2 alignment returning multiple alignments and using longer max insert size limites
# output samtools bam file with only properly aligned paired reads.
bowtie2 -X ${maxInsert} -k ${maxAlignments} --threads ${coreNum} -t -x ${bowtie2GenomeBase} \
    -D 10 -R 2 -L 30 -i S,0,2.50 -N 0 --no-mixed --no-discordant --end-to-end \
    -1 ${tempFolder}/read1_trimmed.fastq.gz -2 ${tempFolder}/read2_trimmed.fastq.gz | \
    samtools view -bh -f 0x003 -o ${tempFolder}/proper_alignment.bam -

# Sort the bam file by name
samtools sort -@ ${coreNum} -m ${memPerCore} -n -o ${tempFolder}/proper_alignment_sortedByName.bam ${tempFolder}/proper_alignment.bam

# Remove reads that fall into off limit regions
bedtools intersect -abam ${tempFolder}/proper_alignment_sortedByName.bam -v -b ${offLimitRegions} > ${tempFolder}/${sampleName}.bam

# Split bam files into 3 bam files
python ${scriptFolder}/split_multiple_alignment.py ${tempFolder}/${sampleName}.bam ${maxAlignments}

# Pull out the uniquely mapped alignments sorted by index and tabulate those
#samtools view ${tempFolder}/proper_alignment_sortedByName.bam | cut -f 1 | uniq -c | awk '{$1=$1;print}' | tr ' ' '\t' > ${tempFolder}/proper_alignment.count
#awk -F '\t' '$1==2 {print $2}' ${tempFolder}/proper_alignment.count > ${tempFolder}/unique_alignment_names.txt
#filterbyname.sh -Xmx16g -eoom in=${tempFolder}/proper_alignment_sortedByName.bam out=${tempFolder}/unique_alignment_sortedByName.bam names=${tempFolder}/unique_alignment_names.txt include=t ow=t

# Pull out the multi-mapped reads and save for later
#awk -F '\t' '$1>2 {print $2}' ${tempFolder}/proper_alignment.count > ${tempFolder}/multiple_alignment_names.txt
#filterbyname.sh -Xmx16g -eoom in=${tempFolder}/proper_alignment_sortedByName.bam out=${tempFolder}/multiple_alignment_sortedByName.bam names=${tempFolder}/multiple_alignment_names.txt include=t ow=t

# Tabulate read count

totalReads=$(( $( zcat ${tempFolder}/read1.fastq.gz | wc -l ) / 4 ))
readsAfterTrim=$(( $( zcat ${tempFolder}/read1_trimmed.fastq.gz | wc -l ) / 4 ))
uniqueReads=$( samtools view ${tempFolder}/${sampleName}.unique_alignments.bam | cut -f1 | uniq | wc -l )
multipleReads1=$( samtools view ${tempFolder}/${sampleName}.multiple_alignments_multiple_genomes.bam | cut -f1 | uniq | wc -l )
multipleReads2=$( samtools view ${tempFolder}/${sampleName}.multiple_alignments_unique_genome.bam | cut -f1 | uniq | wc -l )
multipleReads=$(( multipleReads1 + multipleReads2 ))
echo 'Sample_Name,Total_Fragments,Fragments_After_Trim,Fragments_Aligned_Uniquely,Fragments_Aligned_Multiple_Times' > ${tempFolder}/read_accounting.csv
echo ${sampleName}','${totalReads}','${readsAfterTrim}','${uniqueReads}','${multipleReads} >> ${tempFolder}/read_accounting.csv

# Go through the multimapped reads to assign them using a python script.
date
python ${scriptFolder}/tabulate_alignment_fragment.py -s ${sampleName} ${bowtie2GenomeBase}.fasta ${referenceNameFile} ${tempFolder}/${sampleName}.unique_alignments.bam ${tempFolder}/${sampleName}.multiple_alignments_multiple_genomes.bam ${tempFolder}/${sampleName}.multiple_alignments_unique_genome.bam ${tempFolder}/tabulated_alignment_fragment.csv
date

# Copy all the output files back to S3
aws s3 cp --quiet ${tempFolder}/proper_alignment_sortedByName.bam s3:/${bamOutput}
aws s3 cp --quiet ${tempFolder}/tabulated_alignment_fragment.csv s3:/${relativeAbundanceOutput}
aws s3 cp --quiet ${tempFolder}/read_accounting.csv s3:/${readAccountingOutput}
aws s3 sync --quiet --exclude "*" --include "*.coverage*.csv" ${tempFolder}/ ${s3Root}

pwd
echo "Alignment completed."
ls ${LOCAL}
du -sh ${LOCAL}
rm -rf ${tempFolder}
date
