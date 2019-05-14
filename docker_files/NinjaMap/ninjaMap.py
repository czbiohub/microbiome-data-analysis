#!/usr/bin/env python3

# Each read gets 1 vote for candidate strains. 
# The more strains it maps to the lower the value of it's vote.
# Note: we are not calculating coverage, just read assignments.
# MAJOR ASSUMPTION 1: We will always have a complete genome for organisms (exact strain) that we seek in the sample.

import argparse
import logging
import os
import pysam
import re
import sys
from collections import defaultdict, Counter

usage = """
    USAGE: 
    python ninjaMap.py \
-bam name_sorted_input_bamfile \
-bin tab-delimited file with Col1= contig name and Col2=Bin/Strain name \
-out abundance table output
    """

p = argparse.ArgumentParser(   
    formatter_class=argparse.RawTextHelpFormatter,
    add_help=True,
    usage=argparse.SUPPRESS,
    description="""Description:
This script will calculate the abundance of a strain in a defined microbial community. 
Usage: ninjaMap.py -bam sorted.bam -bin contig_strain_assignments.tsv -out abundance_table_output.tsv
""",
    epilog="""Examples:
python ninjaMap.py -bin contig_names_bin_map.txt -bam Bacteroides-sp-9-1-42FAA/Bacteroides-sp-9-1-42FAA.processed.sortedByCoord.bam -out Bacteroides-sp-9-1-42FAA.sorted.ninjaAbundance.tsv    
""")
p.add_argument('-bam', dest='bamfile', action='store', type=str, required = True,
                help='sorted and indexed bam file.')
p.add_argument('-bin', dest='binmap', action='store', type=str,
                help='tab-delimited file with Col1= contig name and Col2=Bin/Strain name')
p.add_argument('-out', dest='output', action='store', type=str,
                help='abundance table output')
# p.add_argument('-min_pid', dest='min_pid', action='store', type=float, default=0)
# p.add_argument('-min_readq', dest='min_readq', action='store', type=float, default=0)
# p.add_argument('-min_mapq', dest='min_mapq', action='store', type=float, default=10)
# p.add_argument('-min_aln_cov', dest='min_aln_cov', action='store', type=float, default=0)

args = vars(p.parse_args())

# accept bam file (ideally indexed) and contig name -> bin map -> contig_len
# bamfile_name = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/Bacteroides-coprophilus-DSM-18228/Bacteroides-coprophilus-DSM-18228.processed.bam"
# abundance_output_file = 'B_coprophilius.ninjaMap.v1.abundance.tsv'
bamfile_name = args['bamfile']
prefix = os.path.basename(bamfile_name).split('.')[0]

logging.basicConfig(filename=prefix+'.log', 
    filemode='a+', 
    level=logging.DEBUG,
    format='%(asctime)s\t[%(levelname)s]:\t%(message)s')


if not args['binmap']:
    binmap_file = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/contig_names_bin_map.txt"
else:
    binmap_file = args['binmap']

if not args['output']:
    abundance_output_file = prefix +'.ninjaMap_abundance.tsv'
else:
    abundance_output_file = args['output']

logging.info('Started')
logging.info('BAM file: %s', bamfile_name)
logging.info('Bin Map file: %s', binmap_file)
logging.info('Output file: %s', abundance_output_file)

bins = defaultdict()
strains = Counter()
# strain_info = Counter() # strain_info['strain_name'] = genome_size
num_lines = 0
with open(binmap_file, "r") as binmap:
    for line in binmap:
        num_lines += 1
        line=line.rstrip()
        contig_name, strain_name = line.split('\t')
        bins[contig_name] = strain_name
        strains[strain_name] += 1
# print(num_lines)
logging.info('Read %d contigs assigned to %d strains', len(bins.keys()), len(strains.keys()))

## Declaring some default global variables
primary = defaultdict(Counter)
primary_scores = Counter()

escrow = defaultdict(Counter)
escrow_scores = Counter()

read_fractions = defaultdict(Counter)

strain_primary_abundance = Counter()
primary_strain_weights = Counter()

strain_escrow_abundance = Counter()
escrow_strain_weights = Counter()

min_aln_len = 140
num_aln = 0
num_primary_aln = 0
num_escrow_aln = 0

bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
for aln in bamfile.fetch(until_eof=True):
    num_aln += 1
    read = str(aln.qname)
    contig_name = str(aln.reference_name)
    num_primary_aln += 1
    strain_name = bins[contig_name]

    read_fractions['Aligned'][read] += 1
    
    # Also test with PCR duplicates removed as well because it could inflate the read proportion by strain. 
    # Since I'm not looking at coverage, just depth.
    if (dict(aln.tags)['NM'] == 0) and (aln.query_alignment_length == aln.query_length):
        if not aln.is_secondary:
            # true if not primary alignment
            # https://pysam.readthedocs.io/en/latest/api.html?highlight=is_secondary#pysam.AlignedSegment.is_secondary
            # this is a primary alignment
            
            primary[read][strain_name] += 1
            
            if aln.is_duplicate: # primary and duplicate
                read_fractions['Perfect Match Primary Alignment BUT PCR Duplicate'][read] += 1  
            else :
                read_fractions['Perfect Match Primary Alignment'][read] += 1
                # primary_scores['strain_name'] = dict(aln.tags)['AS']
        else:
        # else perfect alignment over at least 140 bases but not primary
            num_escrow_aln += 1
            # The following can also be added, but will make it specific to bowtie2, as the flags used are bowtie2 specific
            # OPTIONAL: only keep secondary if secondary alignment score is within X% of the primary.
                # What was its primary alignment score?
                # What is its secondary alignment score?
                # Difference between primary and secondary alignment scores.
                # Who was the primary alignment with?
                # Who was the secondary alignment with?
            
            escrow[read][strain_name] += 1
            
            if aln.is_duplicate:
                read_fractions['Perfect Match Alignment in Escrow BUT PCR Duplicate'][read] += 1  
            else :
                # escrow_scores['strain_name'] = dict(aln.tags)['AS']
                read_fractions['Perfect Match Alignment in Escrow'][read] += 1
    elif not aln.is_secondary:
        # true if not primary alignment
        # https://pysam.readthedocs.io/en/latest/api.html?highlight=is_secondary#pysam.AlignedSegment.is_secondary
        # this is a primary alignment
        
        escrow[read][strain_name] += 1

        if aln.is_duplicate: # primary and duplicate
            read_fractions['Not Perfect BUT Primary Alignment AND PCR Duplicate'][read] += 1  
        else :
            read_fractions['Not Perfect BUT Primary Alignment'][read] += 1
            # primary_scores['strain_name'] = dict(aln.tags)['AS']
    else:
        escrow[read][strain_name] += 1

        if aln.is_duplicate:
            read_fractions['Not Perfect, Not Primary AND PCR Duplicate'][read] += 1
        else:
            read_fractions['Not Perfect Not Primary'][read] += 1
        # What are these? Write to a new bam file.
#         print(aln)

# print('Total Alignments :' + str(num_aln))
# Actual_Total_Reads = 8356502
# print('Total Reads :'+str(Actual_Total_Reads))

Total_Reads_Aligned = len(read_fractions['Aligned'].keys())
logging.info('Read %d alignments for %d reads', num_aln, Total_Reads_Aligned)

for frac in read_fractions.keys():
    logging.info('%s fraction of all aligned reads : %7.3f%%', frac, len(read_fractions[frac].keys())*100/Total_Reads_Aligned)

strain_primary_abundance = Counter()
primary_strain_weights = Counter()
reads_considered = Counter()
for read in primary.keys():
    # Spread the 1 vote/read evenly amongst all the strains with primary alignment for this read
    ## Should always be 1 strain
    reads_considered[read] += 1
    read_vote_value = 1/len(primary[read].keys())
    for strain in primary[read].keys():
        strain_primary_abundance[strain]+=read_vote_value
    
for strain in strain_primary_abundance.keys():
    primary_strain_weights[strain] = strain_primary_abundance[strain]/Total_Reads_Aligned

# primary_strain_weights        
strain_escrow_abundance = Counter()
escrow_strain_weights = Counter()
for read in escrow.keys():
    if not primary[read]:
        reads_considered[read] += 1
        for strain in escrow[read].keys():
            # Spread 1 vote/read based on the primary weight of strains aligned to
            read_vote_value = 1*primary_strain_weights[strain]
            strain_escrow_abundance[strain] += read_vote_value 
            
for strain in strain_escrow_abundance.keys():
        escrow_strain_weights[strain] = (strain_escrow_abundance[strain]/Total_Reads_Aligned)

# escrow_strain_weights

f = open(abundance_output_file, 'w')
for strain in primary_strain_weights.keys():
    f.write(strain +'\t'+ str(primary_strain_weights[strain] + escrow_strain_weights[strain]) + '\n')

logging.info('%d reads out of %d [~%7.3f%%] aligned reads were used to generate this abundance estimate', 
    len(reads_considered.keys()),
    Total_Reads_Aligned,
    len(reads_considered.keys())*100/Total_Reads_Aligned)

logging.info('Completed')
# Functions for making this more robust
# def bam_is_empty(fn):
#     if os.path.getsize(fn) > 1000000:
#         return False

#     bam = pysam.Samfile(fn, check_sq=False)
#     try:
#         bam.next()
#         return False
#     except StopIteration:
#         return True

# def sort_and_index(file_name, sorted_prefix=None):
#     """ Sorts and indexes a bam file by coordinates.
#     """
#     if sorted_prefix is None:
#         sorted_prefix = file_name.replace('.bam', '') + '_sorted'

#     sorted_name = sorted_prefix + '.bam'
#     pysam.sort(file_name, sorted_prefix)
#     pysam.index(sorted_name)

#     return pysam.Samfile(sorted_name, 'rb')

# def sort_by_name(file_name, sorted_prefix=None):
#     """ Sorts a bam file by the read name, for paired-end
#     """
#     if sorted_prefix is None:
#         sorted_prefix = file_name.replace('.bam', '') + '_namesorted'

#     sorted_name = sorted_prefix + '.bam'
#     pysam.sort('-n', file_name, sorted_prefix)

#     return pysam.Samfile(sorted_name, 'rb')