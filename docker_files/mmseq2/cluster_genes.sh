#!/bin/bash -x

set -euo pipefail

docker container run --rm \
    --volume "$(pwd)/Genes":/app \
    soedinglab/mmseqs2 \
    mmseqs \
        easy-cluster \
        /app/all_complete_genes_combined.fna.gz \
        /app/infantMicrobiome_genes \
        /app/tmpdir \
        --min-seq-id 0.99 \
        -c 0.90 \
        --cov-mode 0 \
        --createdb-mode 0 \
        --dbtype 0 \
        --compressed 1 \
        --cluster-reassign \
        --threads 12 &> mmseq.eclust.log

# # Faster Linclust algorithm
# # Less sensitive
# docker container run --rm \
#     --volume $(pwd)/Genes:/app \
#     soedinglab/mmseqs2 \
#     mmseqs \
#         easy-linclust \
#         /app/all_complete_genes_combined.fna.gz \
#         /app/infantMicrobiome_genes \
#         /app/tmpdir \
#         --min-seq-id 0.99 \
#         -c 0.90 \
#         --cov-mode 0 \
#         --createdb-mode 0 \
#         --dbtype 0 \
#         --threads 12 &> mmseq.elinclust.log
