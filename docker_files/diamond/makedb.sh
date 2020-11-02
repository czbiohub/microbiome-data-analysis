#!/bin/bash -x

set -euo pipefail

docker container run --rm \
    --volume /home/ec2-user/efs/docker/UniRef/UniRef90/20200812:/app/db \
    quay.io/biocontainers/diamond:2.0.2--h56fc30b_0 \
    diamond \
        makedb \
            --in /app/db/uniref90.fasta.gz \
            --db /app/db/uniref90.dmnd


# docker container run --rm \
#     --volume /home/ec2-user/efs/docker/NCBI/20200808:/app/db \
#     quay.io/biocontainers/diamond:2.0.2--h56fc30b_0 \
#     diamond \
#         makedb \
#             --in /app/db/nr.gz \
#             --db /app/db/nr.dmnd \
#             --taxonnodes /app/db/nodes.dmp \
#             --taxonnames /app/db/names.dmp