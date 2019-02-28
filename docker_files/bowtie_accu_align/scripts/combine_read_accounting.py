"""
combine_read_accounting.py

Description: Looks for all files ending with .readAccounting.csv. Each file would
             have been generated from alignment to a combined reference genomes.
             Values represent total number of clusters (paired end reads) 

Author: Brian Yu 2018.09.11
"""

import argparse, os, subprocess

usage = "USAGE: python combine_reed_accounting.py [--s3-dir /czbiohub-brianyu/....] summary_file_location output_filename"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument('--s3-dir', dest='s3_dir', action='store', type=str, default=False)
p.add_argument(dest='summary_file_location', action='store', type=str)
p.add_argument(dest='output_filename', action='store', type=str)

arguments = p.parse_args()

# Read in all the alignment summaries
if arguments.s3_dir:
    # s3 dir, then the folder name should be '/czbiohub-brianyu/......' without / at the end
    subprocess.run(['aws','s3','sync','--exclude','*','--include','*.readAccounting.csv','s3:/'+arguments.s3_dir,arguments.summary_file_location], check=True)

# Now do the common tasks
file_names = os.listdir(arguments.summary_file_location+'/')
file_names = [x for x in file_names if '.readAccounting.csv' in x]
print('Found '+str(len(file_names))+' samples.')

# Read files
header_line = []
reads_line = []
for filename in file_names:
    with open(arguments.summary_file_location+'/'+filename, 'r') as f:
        header_line.append(f.readline()) # first line
        reads_line.append(f.readline()) # second line

# Generate Output
with open(arguments.output_filename, 'w') as output:
    t = output.write(header_line[0])
    t = output.write(''.join(reads_line))
