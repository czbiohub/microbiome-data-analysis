"""
tabulate_alignment_fragment_test.py

reference_list_file: File with two columns used to generate bowtie2 index.
                     Stored with the bowtie2 references and called genome_name_list.csv
unique_alignment_file: bam file with all unique alignments
tabulated_alignment_file: output file including all the genomes and one column with fragments
sample_name: a string that represents the name of the sample

Description: This function is used in accu_align_test.sh. Only the first read alignments are used.
             For unique alignments, the number of bps in inferred fragment length is used.
             This script does not handle multiple alignments anymore.

Author: Brian Yu

Revision History:
2018.09.10 Version 1 had the wrong algorithm. This version uses the unique reads.
2018.09.14 This is the test version to try alignment methods and parameters. 
           This python script is used with the test verion of shell script, where bowtie2 returns
           a single best alignment. This script tabulates the bp fragment size from unique alignments only.
2018.09.15 The next test was to count only reads instead of total fragment lengths. 
           Updated this function to count only reads.
"""

import argparse, pysam
import numpy as np
import pandas as pd

usage = "USAGE: python generate_combined_reference_test.py reference_list_file unique_alignment_file tabulated_alignment_file sample_name"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument(dest='reference_list_file', action='store', type=str)
p.add_argument(dest='unique_alignment_file', action='store', type=str)
# p.add_argument(dest='multiple_alignment_file', action='store', type=str)
p.add_argument(dest='tabulated_alignment_file', action='store', type=str)
p.add_argument(dest='sample_name', action='store', type=str)

arguments = p.parse_args()

# Read in genome names, use the corrected name as the row indices
reference_list = pd.read_csv(arguments.reference_list_file, sep=',', header=0, index_col=1)
reference_list[arguments.sample_name] = 0
coverage = reference_list.drop('file_name', axis=1) # file_name is column header, axis=1 means column

# Process unique alignments
print("Processing Uniquely Aligned Reads")
bamfile = pysam.AlignmentFile(arguments.unique_alignment_file, mode='rb')
# alignment here is a pysam AlignedSegment data structure
counter = 0
for alignment in bamfile.fetch(until_eof=True):
    if alignment.is_paired and alignment.is_read1 and alignment.is_proper_pair:
        ref = bamfile.get_reference_name(alignment.reference_id)
        # coverage.ix[ref.split('_')[0],arguments.sample_name] += float(abs(alignment.template_length))
        coverage.ix[ref.split('_')[0],arguments.sample_name] += 1
        counter += 1
        if counter % 100000 == 0:
            print('.', end='')
            counter = 0
bamfile.close()
print('.')

# For debugging purposes and testing algorithms
coverage.to_csv(arguments.tabulated_alignment_file, index=True, header=True)

# removed the part that handles multiple aligned reads

# Output final tabulated file
coverage.to_csv(arguments.tabulated_alignment_file, index=True, header=True)
