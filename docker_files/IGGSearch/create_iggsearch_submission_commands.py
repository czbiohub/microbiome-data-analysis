"""
seedfile: must be a csv file including all the information for fastq1 fastq2 bamOutput etc.
          the headers must match the required environmental variables in accu_align_v1.sh
output_command: the output file generated with aegea batch commands

Description:    This function generates all the commands needed to submit batch jobs 
                on AWS batch using aegea. Input seedfile required.

Author: Brian Yu

Revision History:
2018.09.10 Created
2018.09.12 Added functionality to check for jobs already done.
2018.09.13 Decrease the number of jobs submitted concurrently and increased wait times.
2019.03.05 Adapted for Sunit's Docker based midas alignment
2019.03.25 Adapted for Sunit's Docker image for igg search alginment, only species
2019.04.01 [SJ] Added retry flag for aegea. 
           [SJ] Simplified script to only produce the aegea commands.
                Use GNU Parallel to submit jobs, like so:
                cat aa_dropout_submission_20190401.sh |\
                    parallel --joblog aa_dropout.joblog --retries 3 --delay 3s \
                    "{}" &> aa_dropout_submission_20190401.log &
2019.04.03 [SJ] Updated finished samples detection
"""

# boto3 is a AWS SDK for python
import argparse, subprocess, os, boto3
import pandas as pd
import re

usage = "USAGE: python create_midas_submission_commands.py [options] seedfile s3output_root output_command"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument(dest='seedfile', action='store', type=str)
p.add_argument(dest='s3output_root', type=str) # include s3://...
p.add_argument(dest='output_command', action='store', type=str)
p.add_argument('-i', '--image', dest='image', action='store', type=str, default='sunitjain/iggsearch:latest')
p.add_argument('-m', '--memory', dest='memory', action='store', type=int, default=32000)
p.add_argument('-c', '--core', dest='vcpus', action='store', type=int, default=8)
p.add_argument('-s', '--storage', dest='storage', action='store', type=int, default=500) # the minimum for AWS is 500
p.add_argument('-q', '--queue', dest='queue', action='store', type=str, default='microbiome-lowPriority')
p.add_argument('-r', '--retry', dest='max_retries', action='store', type=str, default='3')
p.add_argument('--subset', dest='subsetReads', action='store', type=str, default='\'\'')
p.add_argument('--hardTrim', dest='hardTrim', action='store', type=str, default='\'\'')

arguments = p.parse_args()


# Create base command string aegea_batch or aegea_batch_demux
s3_bucket = 's3://czbiohub-microbiome/'
base_string = 'aegea batch submit --retry-attempts '+arguments.max_retries+' --queue '+arguments.queue+' --image '+arguments.image+' --storage /mnt='+str(arguments.storage)+' --memory '+str(arguments.memory)+' --vcpus '+str(arguments.vcpus)+' --command='
command_string1 = '"export coreNum='+str(arguments.vcpus)+'; '
command_string2 = './run_iggsearch.sh"'

# Read in seedfile, column 1 (ie 2) needs to be sampleName
run_samples = pd.read_csv(arguments.seedfile, sep=',', header=0, index_col='sampleName') 
# print(run_samples)

# Read in existing output THIS PART NEEDS TO BE UPDATED
aws_s3_prefix = re.sub(s3_bucket,'',arguments.s3output_root)
# print(aws_s3_prefix)
output_content = boto3.client('s3').list_objects(Bucket=re.sub(r'^s3\:\/\/|\/$','',s3_bucket), Prefix=aws_s3_prefix)
# print(arguments.s3output_root)

finished_samples = []
# boto3's client.list_objects return an 8 length dictionary with 'Content' as the first element if the folder exists
# it returns a length 7 dictionary where 'Content' is not a key
if 'Contents' in output_content.keys():
    account = boto3.client('s3')
    counter = 0
    for sample in run_samples.index:
        # Could add specific Prefix or keys so that we are checking for the existence of a specific file.
        t = account.list_objects(Bucket=re.sub(r'^s3\:\/\/|\/$','',s3_bucket), Prefix=aws_s3_prefix+'/'+sample) # prefix should contain the last slash
        if 'Contents' in t.keys(): # if the sample folder exists
            for i in range(len(t['Contents'])):
                if aws_s3_prefix+'/'+sample+'/job.complete' in t['Contents'][i]['Key']: # if the sample folder contains job.complete
                    # print('Complete')
                    finished_samples.append(sample) # sampleName

# Create command file
with open(arguments.output_command, 'w') as command_file:
    counter = 0
    for sample in run_samples.index:
        # If the sample is not already processed, (subprocess.PIPE must be capitalized)
        if sample not in finished_samples:
            # If you only want a section of the samples analyzed
            sample_specific_file_names = 'export fastq1='+run_samples.loc[sample,'fastq1']+'; export fastq2='+run_samples.loc[sample,'fastq2']+'; export s3OutputPath='+arguments.s3output_root+'/'+sample+'/; export subsetReads='+arguments.subsetReads+'; export hardTrim='+arguments.hardTrim+'; '
            t = command_file.write(base_string+command_string1+sample_specific_file_names+command_string2)
            command_file.write('\n')

