#!/bin/bash

####################################################################
# accu_align_v1.sh
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
#####################################################################

# Prepare the EC2 instance with the ECR image
# /mnt has a volume of defined size mounted
LOCAL=/mnt

# Update utilities repository
cd /utilities
git pull

# Clone microbiome analysis repository
# git clone https://github.com/czbiohub/microbiome-data-analysis.git
cd $LOCAL
ls $LOCAL
adapterFile=$LOCAL/microbiome-data-analysis/misc_files/illumina_adapters.fa
scriptFolder=$LOCAL/microbiome-data-analysis/similarity_aware_alignment/

# Install conda environments
export PATH=~/anaconda/bin:$PATH
# which conda
# which aws
conda config --append channels biobuilds
conda update conda
conda create -q -n alignment-env python=3.6.5 scipy numpy pandas bowtie2 samtools bedtools vcftools pysam sambamba picard fastqc bwa bbmap

# Copy genome reference over
aws s3 sync --quiet s3:/$genomeReferenceDir/ $LOCAL/reference/
bowtie2GenomeBase=$LOCAL/reference/combined_104_reference_genomes
referenceNameFile=$LOCAL/reference/genome_name_list.csv

# Constant definitions for bbduk
trimQuality=25
minLength=50
kmer_value=23
min_kmer_value=11

# Set tempFolder definition
tempFolder=$LOCAL/tmp_$( date +"%Y%m%d_%H%M%S" )
mkdir $tempFolder

# Load environment with all the bioinformatic tools
set +u
source activate alignment-env

echo "Starting to Process Sample: "$sampleName

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet s3:/$fastq1 $tempFolder/read1.fastq.gz
aws s3 cp --quiet s3:/$fastq2 $tempFolder/read2.fastq.gz

# Use bbduk to trim reads, -eoom exits when out of memory
bbduk.sh -Xmx16g -eoom in1=$tempFolder/read1.fastq.gz in2=$tempFolder/read2.fastq.gz out1=$tempFolder/read1_trimmed.fastq.gz out2=$tempFolder/read2_trimmed.fastq.gz ref=$adapterFile ktrim=r k=$kmer_value mink=$min_kmer_value hdist=1 tbo qtrim=rl trimq=$trimQuality minlen=$minLength
	
# bowtie2 alignment returning multiple alignments and using longer max insert size limites
# output samtools bam file with only properly aligned paired reads.
bowtie2 -X $maxInsert -k $maxAlignments --threads $coreNum -t -x $bowtie2GenomeBase -1 $tempFolder/read1_trimmed.fastq.gz -2 $tempFolder/read2_trimmed.fastq.gz | samtools view -bh -f 0x003 -o $tempFolder/proper_alignment.bam -

# Sort the bam file by name
samtools sort -@ $coreNum -m $memPerCore -n -o $tempFolder/proper_alignment_sortedByName.bam $tempFolder/proper_alignment.bam

# Pull out the uniquely mapped alignments sorted by index and tabulate those
samtools view $tempFolder/proper_alignment_sortedByName.bam | cut -f 1 | uniq -c | awk '{$1=$1;print}' | tr ' ' '\t' > $tempFolder/proper_alignment.count
awk -F '\t' '$1==2 {print $2}' $tempFolder/proper_alignment.count > $tempFolder/unique_alignment_names.txt
filterbyname.sh -Xmx16g -eoom in=$tempFolder/proper_alignment_sortedByName.bam out=$tempFolder/unique_alignment_sortedByName.bam names=$tempFolder/unique_alignment_names.txt include=t ow=t

# Pull out the multi-mapped reads and save for later
awk -F '\t' '$1>2 {print $2}' $tempFolder/proper_alignment.count > $tempFolder/multiple_alignment_names.txt
filterbyname.sh -Xmx16g -eoom in=$tempFolder/proper_alignment_sortedByName.bam out=$tempFolder/multiple_alignment_sortedByName.bam names=$tempFolder/multiple_alignment_names.txt include=t ow=t

# Tabulate read count
totalReads=$(( $( zcat $tempFolder/read1.fastq.gz | wc -l ) / 4 ))
readsAfterTrim=$(( $( zcat $tempFolder/read1_trimmed.fastq.gz | wc -l ) / 4 ))
uniqueReads=$( cat $tempFolder/unique_alignment_names.txt | wc -l )
multipleReads=$( cat $tempFolder/multiple_alignment_names.txt | wc -l )
echo 'Sample_Name,Total_Fragments,Fragments_After_Trim,Fragments_Aligned_Uniquely,Fragments_Aligned_Multiple_Times' > $tempFolder/read_accounting.csv
echo $sampleName','$totalReads','$readsAfterTrim','$uniqueReads','$multipleReads >> $tempFolder/read_accounting.csv

# Go through the multimapped reads to assign them using a python script.
python $scriptFolder/tabulate_alignment_fragment_v1.py $referenceNameFile $tempFolder/unique_alignment_sortedByName.bam $tempFolder/multiple_alignment_sortedByName.bam $tempFolder/tabulated_alignment_fragment.csv $sampleName

# Copy all the output files back to S3
aws s3 cp --quiet $tempFolder/proper_alignment_sortedByName.bam s3:/$bamOutput
aws s3 cp --quiet $tempFolder/tabulated_alignment_fragment.csv s3:/$relativeAbundanceOutput
aws s3 cp --quiet $tempFolder/read_accounting.csv s3:/$readAccountingOutput

source deactivate
pwd
echo "Alignment completed."
ls $LOCAL
du -sh $LOCAL
date
