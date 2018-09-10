aegea batch submit --ecr-image aligner --storage /mnt=500 --memory 16000 --vcpus 8 --dry-run \
--command="cd /mnt; \
git clone https://github.com/czbiohub/microbiome-data-analysis.git; \
coreNum=4; \
numPerCore=2G; \
maxInsert=3000; \
maxAlignments=200; \
genomeReferenceDir=/czbiohub-brianyu/Synthetic_Community/Genome_References/Bowtie2Index_090718; \
sampleName=Bacteroides-ovatus-ATCC-8483-20171123-WellE6; \
fastq1=/czbiohub-brianyu/TempFiles/Bacteroides-ovatus-ATCC-8483-20171123-WellE6_S54_R1_001_short.fastq.gz; \
fastq2=/czbiohub-brianyu/TempFiles/Bacteroides-ovatus-ATCC-8483-20171123-WellE6_S54_R2_001_short.fastq.gz; \
bamOutput=/czbiohub-brianyu/TempFiles/Bacteroides-ovatus-ATCC-8483-20171123-WellE6.sortedByName.bam; \
relativeAbundanceOutput=/czbiohub-brianyu/TempFiles/Bacteroides-ovatus-ATCC-8483-20171123-WellE6.alignmentSummary.csv; \
readAccountingOutput=/czbiohub-brianyu/TempFiles/Bacteroides-ovatus-ATCC-8483-20171123-WellE6.readAccounting.csv; \
source /mnt/microbiome-data-analysis/batch_submission_scripts/accu_align_v1.sh"