#!/usr/bin/env python3
"""
General purpose submission script for pipelines that require you to export S3 paths for:
    Output directory    --> 'S3OUTPUTPATH', 
    Input fwd_fastq     --> 'fastq1', and 
    Input rev_fastq     --> 'fastq2'
"""

import s3fs
import argparse
import pandas as pd
import os
import re
import json
import sys

usage = "USAGE: python create_submission_commands.py [options] seedfile s3output_root output_command"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
# Required
p.add_argument('-i','--s3_input', dest='in_s3path', type=str, action='store', required=True)
p.add_argument('-o','--s3_output', dest='out_s3path', type=str, action='store', required=True)
p.add_argument('-a','--commands',dest='aegea_cmd', action='store', type=str, required=True)
p.add_argument('-p','--pipeline', dest='pipeline', action='store', type=str, required=True)
p.add_argument('-n','--project', dest='project', action='store', type=str, required=True)
p.add_argument('-l','--lead', dest='team_lead', action='store', type=str, required=True)
p.add_argument('--config', dest='config_file', action='store', type=str, default = "pipelines.json")

# Optional
p.add_argument('-v','--image_version', dest='image_version', action='store', type=str)
p.add_argument('--image', dest='image_name', action='store', type=str)
p.add_argument('--queue', dest='queue', action='store', type=str)
p.add_argument('--script', dest='script', action='store', type=str)

p.add_argument('-m', '--memory', dest='memory', action='store', type=int, default = 2000)
p.add_argument('-c', '--core', dest='vcpus', action='store', type=int, default = 2)
p.add_argument('-s', '--storage', dest='storage', action='store', type=int, default = 500) # the minimum for AWS is 500
p.add_argument('--retry', dest='max_retries', action='store', type=int, default = 1)

arguments = p.parse_args()
in_s3path = arguments.in_s3path.rstrip("/")
out_s3path = arguments.out_s3path.rstrip("/")
aegea_cmd_file = arguments.aegea_cmd
pipeline = arguments.pipeline.lower()

with open(arguments.config_file) as pipeline_defaults:
    default = json.load(pipeline_defaults)

if default[pipeline]:
    data_mount = default[pipeline]["data_mount"]
    script_loc = default[pipeline]["script"]
    # current setup requires at least 500GB
    storage = max(default[pipeline]["storage"], arguments.storage)
    memory = max(default[pipeline]["memory"], arguments.memory)
    cpu = max(default[pipeline]["cpu"], arguments.vcpus)

    if arguments.queue:
        queue = arguments.queue
    else:
        queue = default[pipeline]["queue"]

    if queue == "microbiome-highPriority":
        retries = min(1, arguments.max_retries)
    else:
        retries = min(3, arguments.max_retries)

    if arguments.image_name:
        image = arguments.image_name
    else:
        image = default[pipeline]["image"]

    if arguments.image_version:
        version = arguments.image_version
    else:
        version = default[pipeline]["version"]
else:
    print(f"This script has not been configured to be used for {pipeline}. Please update {arguments.config_file} and try again.")
    sys.exit(1)


fs = s3fs.S3FileSystem(anon=False)

def submit_job (name, fwd, rev, s3output = out_s3path,
                job_queue = queue, img_version = version, docker_image = image,
                job_storage = storage, job_cpu = cpu, job_memory = memory,
                job_data_mount = data_mount, job_script = script_loc):

    aegea_cmd = ''
    s3output_path = f"{s3output}/{name}"
    complete = fs.exists(f'{s3output_path}/job.complete')
    if not complete:
        execute_cmd = f'export coreNum={job_cpu};export fastq1={fwd};export fastq2={rev};export S3OUTPUTPATH={s3output_path}; {job_script}' 
        aegea_cmd = f'aegea batch submit --retry-attempts 1 --queue {job_queue} --image {docker_image}:{img_version} --storage {job_data_mount}={job_storage} --memory {job_memory} --vcpus {job_cpu} --command=\'{execute_cmd}\''
    
    return(aegea_cmd)

df = pd.DataFrame()

df["FilePaths"] = fs.glob(in_s3path + '/*fastq.gz')
df["sampleName"] = df["FilePaths"].str.extract(r'([a-zA-Z0-9_\-\.]+)_S\d+_R[12]_\d+\.fastq\.gz')
df["Orientation"] = df["FilePaths"].str.extract(r'[a-zA-Z0-9_\-\.]+_S\d+_(R[12])_\d+\.fastq\.gz')
df["FilePaths"] = df["FilePaths"].apply(lambda x: 's3://' + x)

seed_df = df.pivot(index='sampleName',columns='Orientation',values='FilePaths')
seed_df.columns.name = None
seed_df = seed_df.reset_index()
# seed_df ["commands"] = seed_df.apply(lambda row: submit_job(name = row['sampleName'], fwd = row['R1'], rev = row['R2']), axis=1)

commands = seed_df.apply(lambda row: submit_job(name = row['sampleName'], fwd = row['R1'], rev = row['R2']), axis=1)
commands.to_csv(aegea_cmd_file, header = False, index = False)