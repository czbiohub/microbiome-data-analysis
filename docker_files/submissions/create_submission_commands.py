#!/usr/bin/env python3
"""
General purpose submission script for pipelines that require you to export S3 paths for:
    Output directory    --> 'S3OUTPUTPATH', 
    Input fwd_fastq     --> 'R1', and 
    Input rev_fastq     --> 'R2'
To use this script with more pipleines, update the pipelines.json file.
"""

import argparse
import json
import os
import re
import sys

import pandas as pd
import s3fs

usage = "python create_submission_commands.py   --seed seedfile.txt \
                                                --s3_output s3://czbiohub-microbiome/My_Lab/My_Name/My_Project/ \
                                                --commands aegea_launch_commands.sh \
                                                --pipeline midas \
                                                --config pipeline.json"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
# Required
p.add_argument('-o','--s3_output', dest='out_s3path', action='store', type=str, required=True,
                help = 'S3 Path to the output directory to store the results')
p.add_argument('-a','--commands',dest='aegea_cmd', action='store', type=str, required=True,
                help = 'name of the file where the aegea commands for each input will be stored')
p.add_argument('-p','--pipeline', dest='pipeline', action='store', type=str, required=True,
                help = 'name of the pipeline to execute. Should be a valid pipeline, present in the pipelines.json file.')
p.add_argument('-c','--config', dest='config_file', action='store', type=str, default = "pipelines.json", required=True,
                help = 'path to the pipelines.json file')

#Pick one
p.add_argument('-s','--seed', dest='seedfile', action='store', type=str,
                help = '(preferred) comma-delimited file containing 3 columns, named: "sampleName","R1" and "R2"')
p.add_argument('-i','--s3_input', dest='in_s3path', action='store', type=str,
                help = 'S3 Path to paired fastq files to be used as inputs in the pipeline')

# Optional
## AWS Logistics
p.add_argument('--to_aws', dest='xfer', action='store_true')
p.add_argument('--yes', dest='no_confirmation', action='store_true')
p.add_argument('--image', dest='image_name', action='store', type=str)
p.add_argument('--image_version', dest='image_version', action='store', type=str)
p.add_argument('--queue', dest='queue', action='store', type=str)
p.add_argument('--script', dest='script', action='store', type=str)
p.add_argument('--extra', dest='extraList', action='append', type=str, default = [],
                help = 'Use this flag to add custom export variables. Can be used multiple times in a single command to add multiple variables. You will be asked for a confirmation (unless --yes specified).')
p.add_argument('--force', dest='overwrite', action='store_true', default=False,
                help = 'Overwrite existing results')

## Hardware
p.add_argument('--memory', dest='memory', action='store', type=int, default = 2000)
p.add_argument('--core', dest='vcpus', action='store', type=int, default = 2)
p.add_argument('--storage', dest='storage', action='store', type=int, default = 500) # the minimum for AWS is 500
p.add_argument('--retry', dest='max_retries', action='store', type=int)

## Tagging
p.add_argument('-n','--project', dest='project', action='store', type=str,
                help = 'name of the project that this set of jobs will belong to')
p.add_argument('-l','--lead', dest='team_lead', action='store', type=str,
                help = 'name of the team lead / PI responsible for this set of jobs')

arguments = p.parse_args()
out_s3path = arguments.out_s3path.rstrip("/")
aegea_cmd_file = arguments.aegea_cmd
pipeline = arguments.pipeline.lower()

with open(arguments.config_file) as pipeline_defaults:
    default = json.load(pipeline_defaults)

# Set Defaults
if default[pipeline]:
    data_mount = default[pipeline]["data_mount"]
    script_loc = default[pipeline]["script"]
    # current setup requires at least 500GB
    storage = default[pipeline]["storage"]
    memory = default[pipeline]["memory"]
    cpu = default[pipeline]["cpu"]
    queue = default[pipeline]["queue"]
    image = default[pipeline]["image"]
    version = default[pipeline]["version"]
else:
    print(f"[FATAL] This script has not been configured to be used for {pipeline}.\nPlease update {arguments.config_file} and try again.")
    sys.exit(1)

# Overriding Defaults
if arguments.queue:
    queue = arguments.queue

if arguments.max_retries:
    retries = arguments.max_retries
elif default[pipeline]["retries"]:
    retries = default[pipeline]["retries"]
else:
    retries = 1

if arguments.image_name:
    image = arguments.image_name

if arguments.image_version:
    version = arguments.image_version

if arguments.storage:
    storage = max(storage,arguments.storage)

if arguments.memory:
    memory = max(memory, arguments.memory)

if arguments.vcpus:
    cpu = max(cpu, arguments.vcpus)

overwrite = False
if arguments.overwrite:
    overwrite = True

extra = None
if len(arguments.extraList) > 0:
    extraList = map(lambda x: str(x).rstrip(r';$'), arguments.extraList)
    extra = ''.join([f'export {e};' for e in extraList])
    if not arguments.no_confirmation:
        while True:
            confirm = input(f'''
        The following statement(s) will be added, as is to the export variables:

        {extra}

        Are you sure you wish to continue? [Y/n]
        ''')
            if confirm == 'n':
                sys.exit(0)
            elif confirm == 'Y':
                print('Confirmed. Processing inputs ...')
                break
            else:
                print(f'\n[Invalid Response] "{confirm}" is not a valid response. Please select either "Y" or "n" (case-sensitive)\nTry again:')
else:
    extra = ''

fs = s3fs.S3FileSystem(anon=False)

def submit_job (name, fwd, rev, s3output = out_s3path,
                job_queue = queue, img_version = version, docker_image = image,
                job_storage = storage, job_cpu = cpu, job_memory = memory, job_attempts = retries,
                job_data_mount = data_mount, job_script = script_loc, pipeline = pipeline, extra = extra, overwrite = overwrite):

    aegea_cmd = None
    memoryPerCore = f"{int(job_memory/(job_cpu * 1000))}G"
    s3output_path = f"{s3output}/{name}"
    if overwrite:
        complete = False
    else:
        complete = fs.exists(f'{s3output_path}/job.complete')

    if not complete:
        execute_cmd = f'{extra}export coreNum={job_cpu};export memPerCore={memoryPerCore};export fastq1={fwd};export fastq2={rev};export S3OUTPUTPATH={s3output_path}; {job_script}' 
        aegea_cmd = f'aegea batch submit --retry-attempts {job_attempts} --name {pipeline}__{name} --queue {job_queue} --image {docker_image}:{img_version} --storage {job_data_mount}={job_storage} --memory {job_memory} --vcpus {job_cpu} --command=\'{execute_cmd}\''

    return aegea_cmd

seed_df = pd.DataFrame()
if arguments.seedfile:
    try:
        # Comma sep
        seed_df = pd.read_csv(arguments.seedfile)
        seed_df["sampleName"]
    except KeyError as e:
        # Tab sep
        seed_df = pd.read_table(arguments.seedfile)
        seed_df["sampleName"]
    except KeyError as e:
        print(f"[FATAL] {e}")
        print(df.sample(10))
        sys.exit(1)
    finally:
        print(f'Sample:\n{seed_df.head(5)}')
elif arguments.in_s3path:
    in_s3path = arguments.in_s3path.rstrip("/")
    df = pd.DataFrame()
    df["FilePaths"] = fs.glob(in_s3path + '/*fastq.gz')
    
    try:
        df["sampleName"] = df["FilePaths"].str.extract(r'([a-zA-Z0-9_\-\.]+)_S\d+_L\d+_R[12]_\d+\.fastq\.gz')
        df["Orientation"] = df["FilePaths"].str.extract(r'[a-zA-Z0-9_\-\.]+_S\d+_L\d+_(R[12])_\d+\.fastq\.gz')
        df["FilePaths"] = df["FilePaths"].apply(lambda x: 's3://' + x)
        seed_df = df.pivot(index='sampleName',columns='Orientation',values='FilePaths')
    except ValueError as e:
        df["sampleName"] = df["FilePaths"].str.extract(r'([a-zA-Z0-9_\-\.]+)_S\d+_R[12]_\d+\.fastq\.gz')
        df["Orientation"] = df["FilePaths"].str.extract(r'[a-zA-Z0-9_\-\.]+_S\d+_(R[12])_\d+\.fastq\.gz')
        # df["FilePaths"] = df["FilePaths"].apply(lambda x: 's3://' + x)
        seed_df = df.pivot(index='sampleName',columns='Orientation',values='FilePaths')
    except e:
        print(f"[FATAL] {e}")
        print(df.sample(10))
        sys.exit(1)
    finally:
        print(f'Sample:\n{seed_df.head(5)}')

    seed_df.columns.name = None
    seed_df = seed_df.reset_index()
    seedfile = f'{os.path.splitext(aegea_cmd_file)[0]}.seedfile.csv'
    seed_df.to_csv(seedfile, index = False)
    # seed_df ["commands"] = seed_df.apply(lambda row: submit_job(name = row['sampleName'], fwd = row['R1'], rev = row['R2']), axis=1)
else:
    print('[FATAL] At least one form of input required. Please either submit a seedfile or a s3 folder with paired fastq.gz files')
    sys.exit(1)

# Handle special characters in Sample Names
seed_df["sampleName"] = seed_df["sampleName"].map(lambda x: re.sub(r'\W+', '-', x))

commands = seed_df.apply(lambda row: submit_job(name = row['sampleName'], fwd = row['R1'], rev = row['R2']), axis=1).dropna(axis='index')
commands.to_csv(aegea_cmd_file, header = False, index = False)

if arguments.xfer and (commands.shape[0] > 0):
    aegea_cmd_filename = os.path.basename(aegea_cmd_file)
    fs.put(aegea_cmd_file, f'{out_s3path}/{aegea_cmd_filename}')
