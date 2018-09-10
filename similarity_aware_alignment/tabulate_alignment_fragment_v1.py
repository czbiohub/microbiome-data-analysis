"""
tabulate_alignment_fragment.py

reference_list_file: File with two columns used to generate bowtie2 index.
                     Stored with the bowtie2 references and called genome_name_list.csv
unique_alignment_file: bam file with all unique alignments
multiple_alignment_file: bam file with all reads that mapped multiple times
tabulated_alignment_file: output file including all the genomes and one column with fragments
sample_name: a string that represents the name of the sample

Description: This function is used in accu_align.sh. Only the first read alignments are used.
             For unique alignments, the number of bps in inferred fragment length is used.
             For multiple alignments, the number of bps are proportionally assigned.
"""

import argparse, pysam
import pandas as pd

usage = "USAGE: python generate_combined_reference_v1.py reference_list_file unique_alignment_file multiple_alignment_file tabulated_alignment_file sample_name"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument(dest='reference_list_file', action='store', type=str)
p.add_argument(dest='unique_alignment_file', action='store', type=str)
p.add_argument(dest='multiple_alignment_file', action='store', type=str)
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
        coverage.ix[ref.split('_')[0],arguments.sample_name] += float(abs(alignment.template_length))
        counter += 1
        if counter % 100000 == 0:
            print('.', end='')
            counter = 0
bamfile.close()
print('.')

# For debugging purposes
# coverage.to_csv(arguments.tabulated_alignment_file, index=True, header=True)

# Process multiple alignments
print("Processing Reads With Multiple Alignments")
# Create another df with same index
tempCov = pd.DataFrame(index=coverage.index)
tempCov[arguments.sample_name] = 0
current_read_name = ''
bamfile = pysam.AlignmentFile(arguments.multiple_alignment_file, mode='rb')
for alignment in bamfile.fetch(until_eof=True):
    if alignment.is_paired and alignment.is_read1 and alignment.is_proper_pair:
        # If it's a new read
        if alignment.query_name != current_read_name:
            # Split the read and add to the original coverage df
            if tempCov[arguments.sample_name].sum() > 0:
                coverage[arguments.sample_name] = coverage[arguments.sample_name] + tempCov[arguments.sample_name] / tempCov[arguments.sample_name].sum()
            # Update the variables to handle the next read
            tempCov[arguments.sample_name] = 0
            current_read_name = alignment.query_name
        ref = bamfile.get_reference_name(alignment.reference_id)
        tempCov.ix[ref.split('_')[0],arguments.sample_name] += float(abs(alignment.template_length))
        counter += 1
        if counter % 100000 == 0:
            print('.', end='')
            counter = 0
bamfile.close()
print('.')

# Output final tabulated file
coverage.to_csv(arguments.tabulated_alignment_file, index=True, header=True)
