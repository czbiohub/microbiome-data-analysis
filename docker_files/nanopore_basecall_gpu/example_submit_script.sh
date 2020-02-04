# Potentially make a queue that runs on p3.2xlarge instances with 
# the Deep Learning AMI (v24) pre-installed (https://aws.amazon.com/marketplace/pp/B077GCH38C, ami-0ddba16a97b1dcda5

aegea batch submit --queue no_queue_made_yet --image bmerrill9/nanopore_basecall_demux_gpu:latest --storage /mnt=500 --memory 60000 --vcpus 8 --command="\
export coreNum=8; \
export mem_mb=60000; \
export NANOPORE_RUNPATH=s3://czb-seqbot/nanopore; \
export RUN_NAME=181219_MN19452_0007_FAK54010-prevotella-copri-2379; \
export barcoded=no; \
export flowcell=FLO-MIN106; \
export kit=SQK-LSK109; \
export barcodekit=EXP-NBD104; \
export DATA_DEST=s3://czbiohub-microbiome/Sonnenburg_Lab/Strain_Genome_Sequencing/Nanopore_Basecall_Demux; \
./nanopore_basecall_demux.sh"
