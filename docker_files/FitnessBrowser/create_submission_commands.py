#!/usr/bin/env python3

import s3fs
import argparse
import pandas as pd
import os

fs = s3fs.S3FileSystem(anon=False)

# in_s3path
# out_s3path
# sampleSheet
# job_submissions
usage = "USAGE: python create_submission_commands.py -h"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument('--sample_sheet', dest='sampleSheet', action='store', type=str, required=True)
p.add_argument('--s3_input_dir',dest='in_s3path', type=str, action='store', required=True)
p.add_argument('--s3_output_dir',dest='out_s3path', action='store', type=str, required=True)
p.add_argument('--commands',dest='job_submissions', action='store', type=str, required=True)
p.add_argument('--s3_job_dir',dest='job_s3path', type=str, action='store', required=False, 
    default ="s3://czbiohub-microbiome/Job_Submissions/BarSeq/")
p.add_argument('-i', '--image', dest='image', action='store', type=str, default='sunitjain/fibo:latest')
p.add_argument('-m', '--memory', dest='memory', action='store', type=int, default=8000)
p.add_argument('-c', '--core', dest='vcpus', action='store', type=int, default=4)
p.add_argument('-s', '--storage', dest='storage', action='store', type=int, default=500) # the minimum for AWS is 500
p.add_argument('-q', '--queue', dest='queue', action='store', type=str, default='microbiome-highPriority')
p.add_argument('-r', '--retry', dest='max_retries', action='store', type=str, default='3')

arguments = p.parse_args()
in_s3path = arguments.in_s3path.rstrip("/")
out_s3path = arguments.out_s3path.rstrip("/")
sampleSheet = arguments.sampleSheet
job_submissions = arguments.job_submissions
queue = arguments.queue
memory = arguments.memory
cpu = arguments.vcpus
storage = arguments.storage

version = "latest"
data_mount ="/data"
script_loc = "/mnt/run_multiCodes.sh"

def submit_job (s3input,s3output,index_name,
                job_queue = queue, img_version = version,
                job_storage = storage, job_cpu = cpu, job_memory = memory,
                job_data_mount = data_mount, job_script = script_loc):
    
    execute_cmd = f'export coreNum={job_cpu};export S3INPUTPATH={s3input};export S3OUTPUTPATH={s3output};export INDEX_NAME={index_name}; {job_script}' 
    aegea_cmd = f'aegea batch submit --retry-attempts 1 --queue {job_queue} --image sunitjain/fibo:{img_version} --storage {job_data_mount}={job_storage} --memory {job_memory} --vcpus {job_cpu} --command=\'{execute_cmd}\''
    
    return(aegea_cmd)

df = pd.read_csv(sampleSheet, skiprows = 20,  usecols = ["Sample_Name","Index_ID"])
df["S3Output"] = out_s3path + '/' + df.Sample_Name
df["Index_Name"] = df.Index_ID.str.split("_", expand = True)[2]
df["S3Input"] = df["Sample_Name"].apply(lambda x: 's3://' + fs.glob(in_s3path + '/'+ x +'*.fastq.gz')[0])
df["finished"] = df["Sample_Name"].apply(lambda x: fs.exists(out_s3path + '/'+ x +'/job.complete'))

jobs_remaining = df[df.finished == False]

# Files for BarSeqTest.pl
index_df = jobs_remaining[["Sample_Name","Index_Name","S3Input"]]
index_df.to_csv(f'{os.path.basename(job_submissions)}.index.csv', header = False, index = False)

#Copy to S3

# Create aegea submission commands
commands = jobs_remaining.apply(lambda row: submit_job(row['S3Input'], row['S3Output'],row['Index_Name']), axis=1)
commands.to_csv(job_submissions, header = False, index = False)