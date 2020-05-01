import argparse
import json
import os
import re
import sys

import pandas as pd
import s3fs

seed_file = sys.argv[1]
aegea_cmd_file = sys.argv[2]

# Read the seed_file
df = pd.read_csv(seed_file)
df.head()

# Set pipeline default values
job_count = 0
job_cpu = 16
job_memory = 64000
memoryPerCore = int(job_memory/((job_cpu + 1)*1000))
job_attempts = 1
job_data_mount = '/mnt'
job_storage = 500
job_queue = 'microbiome-highPriority'

pipeline = 'bowtie2'
job_script = './run_bowtie2.sh'
docker_image = 'sunitjain/bowtie2'
img_version = '20200212111032'
skip_qc = 'false'

with open(aegea_cmd_file, 'w') as fh:
    for rowidx, row in df.iterrows():
        s3dbfasta = row['S3DBFASTA']
        fwd = row['fastq1']
        rev = row['fastq2']
        s3output_path = row['S3OUTPUTPATH']
        job_count += 1
        name = os.path.splitext(os.path.basename(s3dbfasta))[0]

        execute_cmd = f'export SKIP_QC={skip_qc};export coreNum={job_cpu};export memPerCore={memoryPerCore};export S3DBFASTA={s3dbfasta};export fastq1={fwd};export fastq2={rev};export S3OUTPUTPATH={s3output_path}; {job_script}' 
        aegea_cmd = f'aegea batch submit --retry-attempts {job_attempts} --name {job_count}__{pipeline}__{name} --queue {job_queue} --image {docker_image}:{img_version} --storage {job_data_mount}={job_storage} --memory {job_memory} --vcpus {job_cpu} --command=\'{execute_cmd}\''

        fh.write(f'{aegea_cmd}\n')
