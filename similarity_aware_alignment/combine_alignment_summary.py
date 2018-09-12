"""
combine_alignment_summary.py

Description: Looks for all files ending with .alignmentSummary.csv. Each file would
             have been generated from alignment to a combined reference genomes.
             Values represent total number of bps from sequenced fragments that align
             to a particular genome. Normalize by the genome size. Then concatenate.

Author: Brian Yu 2018.09.10
"""

import argparse, os, subprocess
import pandas as pd

usage = "USAGE: python combine_alignment_summary.py [--s3-dir /czbiohub-brianyu/....] reference_stats_file summary_file_location output_filename"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument('--s3-dir', dest='s3_dir', action='store', type=str, default=False)
p.add_argument(dest='reference_stats_file', action='store', type=str)
p.add_argument(dest='summary_file_location', action='store', type=str)
p.add_argument(dest='output_filename', action='store', type=str)

arguments = p.parse_args()

# Read in genome names, use the corrected name as the row indices
# The column that contains genome size is size_of_genome and it's in bps
reference_stats = pd.read_csv(arguments.reference_stats_file, sep=',', header=0, index_col=1)

# Read in all the alignment summaries
if arguments.s3_dir:
    # s3 dir, then the folder name should be '/czbiohub-brianyu/......' without / at the end
    subprocess.run(['aws','s3','sync','--exclude','*','--include','*.alignmentSummary.csv','s3:/'+arguments.s3_dir,arguments.summary_file_location], check=True)

# Now do the common tasks
summary_list = []
file_names = os.listdir(arguments.summary_file_location+'/')
file_names = [x for x in file_names if '.alignmentSummary.csv' in x]
print('Found '+str(len(file_names))+' samples.')
for filename in file_names:
    temp_df = pd.read_csv(arguments.summary_file_location+'/'+filename, sep=',', header=0, index_col=0)
    sampleName = list(temp_df)[0] # always the first element
    k = temp_df.loc[temp_df.index, sampleName] / reference_stats.loc[reference_stats.index, 'size_of_genome']
    summary_list.append(k.to_frame(name=sampleName)) # k is actually a pandas series
    # print(k.reset_index(name=sampleName))

# Concat dataframes
pd.concat(summary_list, axis=1, join='outer').to_csv(arguments.output_filename, index=True, header=True)
