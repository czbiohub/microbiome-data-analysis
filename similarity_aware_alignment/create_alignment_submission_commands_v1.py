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
"""

import argparse, subprocess, os
import pandas as pd

usage = "USAGE: python create_alignment_submission_commands_v1.py seedfile output_command"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument(dest='seedfile', action='store', type=str)
p.add_argument(dest='output_command', action='store', type=str)

arguments = p.parse_args()

# Create base command string
base_string = 'aegea batch submit --ecr-image aligner --storage /mnt=500 --memory 32000 --vcpus 8 --command='
command_string1 = '"cd /mnt; git clone https://github.com/czbiohub/microbiome-data-analysis.git; coreNum=8; memPerCore=2G; maxInsert=3000; maxAlignments=200; genomeReferenceDir=/czbiohub-brianyu/Synthetic_Community/Genome_References/Bowtie2Index_090718; '
command_string2 = 'source /mnt/microbiome-data-analysis/similarity_aware_alignment/accu_align_v1.sh"'
# command_string2 = 'source /mnt/microbiome-data-analysis/similarity_aware_alignment/accu_align_test.sh"'

# Read in seedfile, column 1 (ie 2) needs to be sampleName
run_samples = pd.read_csv(arguments.seedfile, sep=',', header=0, index_col='sampleName') 
# print(run_samples)

# Read in existing output
tmp_file = 'tmp_aws_files.txt' # in the current folder, but later deleted
aws_s3_output_folder = run_samples.loc[run_samples.index[0],'relativeAbundanceOutput'].split('/')
aws_s3_output_folder = 's3:/' + '/'.join(aws_s3_output_folder[0:-1]) + '/'
print(aws_s3_output_folder)
with open(tmp_file, 'w') as f:
    t = subprocess.call(['aws','s3','ls',aws_s3_output_folder], stdout=f)
finished_samples = []
if t==0 and os.path.exists(tmp_file) and os.path.getsize(tmp_file) > 0:
    with open(tmp_file, 'r') as f:
        for l in f:
            if 'alignmentSummary.csv' in l:
                t = l.split(' ')[-1] # last element
                finished_samples.append(t.split('.')[0]) # sampleName
else:
    t = subprocess.call(['rm',tmp_file])
# print(finished_samples)

# Create command file
with open(arguments.output_command, 'w') as command_file:
    counter = 0
    for sample in run_samples.index:
        # If the sample is not already processed, (subprocess.PIPE must be capitalized)
        if sample not in finished_samples:
            # If you only want a section of the samples analyzed
            # if counter >= 600:
            sample_specific_file_names = 'fastq1='+run_samples.loc[sample,'fastq1']+'; fastq2='+run_samples.loc[sample,'fastq2']+'; bamOutput='+run_samples.loc[sample,'bamOutput']+'; relativeAbundanceOutput='+run_samples.loc[sample,'relativeAbundanceOutput']+'; readAccountingOutput='+run_samples.loc[sample,'readAccountingOutput']+'; '
            t = command_file.write('echo "Submitting Sample '+sample+'"\n')
            t = command_file.write(base_string+command_string1+'sampleName='+sample+'; '+sample_specific_file_names+command_string2)
            if counter % 8 == 7:
                t = command_file.write('\n\nsleep 600\n\n')
            else:
                t = command_file.write('\n\nsleep 150\n\n')
            counter += 1

# Remove temp file
if os.path.exists(tmp_file):
    t = subprocess.call(['rm',tmp_file])

