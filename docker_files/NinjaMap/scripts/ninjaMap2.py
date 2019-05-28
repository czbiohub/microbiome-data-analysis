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

###############################################################################
# Functions for making this more robust
###############################################################################
def bam_is_empty(fn):
    if os.path.getsize(fn) > 1000000:
        return False

    bam = pysam.Samfile(fn, check_sq=False)
    try:
        bam.next()
        return False
    except StopIteration:
        return True

def sort_and_index(file_name, sorted_prefix=None, cores=1):
    """ Sorts and indexes a bam file by coordinates.
    """
    if sorted_prefix is None:
        sorted_prefix = file_name.replace('.bam', '') + '_sorted'

    sorted_name = sorted_prefix + '.bam'
    pysam.sort('-@',cores, file_name, sorted_prefix)
    pysam.index(sorted_name)

    return pysam.Samfile(sorted_name, 'rb')

def sort_by_name(file_name, sorted_name=None, cores=1):
    """ Sorts a bam file by the read name, for paired-end
    """
    if sorted_name is None:
        sorted_name = file_name.replace('.bam', '') + '_namesorted.bam'

    pysam.sort('-@',cores,'-n', '-o' , sorted_name, file_name)

    return pysam.Samfile(sorted_name, 'rb')

###############################################################################
###############################################################################

usage = """
    USAGE: 
    python ninjaMap.py \
-bam input_bamfile \
-bin tab-delimited file with Col1= contig name and Col2=Bin/Strain name \
-out abundance table output
    """

p = argparse.ArgumentParser(   
    formatter_class=argparse.RawTextHelpFormatter,
    add_help=True,
    usage=argparse.SUPPRESS,
    description="""Description:
This script will calculate the abundance of a strain in a defined microbial community. 
Usage: ninjaMap.py -bam name_sorted.bam -bin contig_strain_assignments.tsv -out abundance_table_output.tsv
""",
    epilog="""Examples:
python ninjaMap.py -bin contig_names_bin_map.txt -bam Bacteroides-sp-9-1-42FAA/Bacteroides-sp-9-1-42FAA.processed.sortedByCoord.bam -out Bacteroides-sp-9-1-42FAA.sorted.ninjaAbundance.tsv    
""")
p.add_argument('-bam', dest='bamfile', action='store', type=str, required = True,
                help='name sorted bam file.')
p.add_argument('-bin', dest='binmap', action='store', type=str,
                help='tab-delimited file with Col1= contig name and Col2=Bin/Strain name')
p.add_argument('-out', dest='output', action='store', type=str,
                help='abundance table output')
p.add_argument('-log', dest='logfile', action='store', type=str,
                help='job log and run summary')
p.add_argument('-threads', dest='threads', action='store', type=str, default=1,
                help='number of threads available for this job and subprocesses')
args = vars(p.parse_args())

# DEBUG
# bamfile_name = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/Bacteroides-coprophilus-DSM-18228/Bacteroides-coprophilus-DSM-18228.processed.bam"
# abundance_output_file = 'B_coprophilius.ninjaMap.v1.abundance.tsv'
bamfile_name = args['bamfile']
prefix = os.path.basename(bamfile_name).split('.')[0]

if not args['logfile']:
    logfile = prefix +'.ninjaMap.log'
else:
    logfile = args['logfile']

logging.basicConfig(filename=logfile, 
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

###############################################################################
# Classes
###############################################################################
class Read:
    # self.is_primary = not aln_obj.is_secondary
    # determine if read is primary or escrow?
    # calculate the strain distribution based on the primary votes.
    # calculate the escrow vote based on the primary vote strain distributions.
    total_alignments_considered = 0
    num_read_alignments = 0
    num_perfect_alignments = 0
    num_escrow = 0
    best_aln_length = 0
    best_aln_score = 0
    escrow = defaultdict(Counter())
    perfect_alignment = defaultdict(Counter())
    abundance = defaultdict(Counter())

    def __init__(self, aln_obj):
        if aln_obj == None:
            self.name = ''
        else:
            self.is_proper_pair = aln_obj.is_proper_pair()
            self.name = aln_obj.qname
            self.ref_name = aln_obj.reference_name
            self.query_aln_length = aln_obj.query_alignment_length
            self.query_length = aln_obj.query_length
            self.aln_cov = align_len/float(self.query_len)
            self.quality = aln_obj.mapping_quality
            self.alignment_score = dict(aln_obj.tags)['AS']
            self.mate_score = dict(aln_obj.tags)['YS']
            self.mismatches = dict(aln_obj.tags)['XM']
            self.gap_open = dict(aln_obj.tags)['XO']
            self.gap_ext = dict(aln_obj.tags)['XG']
            self.is_dup = aln_obj.is_duplicate
            self.is_supp = aln_obj.is_supplementary
            self.is_primary = not(aln_obj.is_secondary)
            self.edit_dist = dict(aln_obj.tags)['NM']
            self.perc_id = 100*(self.align_len-self.edit_dist)/float(self.align_len)
                       
            self.strain = bins[self.ref_name]
            self.num_read_alignments += 1

            if self.query_aln_length > Read.best_aln_length:
                Read.best_aln_length = self.query_aln_length

            if self.alignment_score > self.best_aln_score:
                self.best_aln_score = self.alignment_score

    def keep_alignment(self):
        if (self.edit_dist == 0) and (self.query_aln_length == self.query_length):
            Read.total_alignments_considered += 1
            self.perfect_alignment[self.strain] += 1
            return True
        else:
            return False
    
    def to_escrow(self):
        self.num_strains = len(self.perfect_alignment.keys())

        if self.num_strains == 1:
            return False
        else:
            return True
        
    def tally_votes(self):
        if self.to_escrow():
            self.add_to_escrow()
            return 0
        else:
            # True primary; consider a full vote.
            self.add_to_primary()
            return 1

    def add_to_escrow(self):
        Read.num_escrow += 1
        for strains in self.perfect_alignment.keys():
            Read.escrow[self.name][strains] += 1

    def calculate_escrow_votes(self):
        # The following can also be added, but will make it specific to bowtie2, as the flags used are bowtie2 specific
            # OPTIONAL: only keep secondary if secondary alignment score is within X% of the primary.
                # What was its primary alignment score?
                # What is its secondary alignment score?
                # Difference between primary and secondary alignment scores.
                # Who was the primary alignment with?
                # Who was the secondary alignment with?
        pass

    def calculate_final_weights(self, vote_value):
        Read.abundance[self.strain] += vote_value

    # possible class method
    def add_to_primary(self):
        Read.abundance[self.strain] += 1
            
###############################################################################
###############################################################################

logging.info('Started')
logging.info('BAM file: %s', bamfile_name)
logging.info('Bin Map file: %s', binmap_file)
logging.info('Output file: %s', abundance_output_file)

strains = Counter()
bins = defaultdict()
num_strains = 0
abundance = defaultdict(Counter())
# strain_info = Counter() # strain_info['strain_name'] = genome_size
num_lines = 0
with open(binmap_file, "r") as binmap:
    for line in binmap:
        num_lines += 1
        line=line.rstrip()
        contig_name, strain_name = line.split('\t')
        bins[contig_name] = strain_name
        num_strains += 1
# print(num_lines)
logging.info('Read %d contigs assigned to %d strains', len(bins.keys()), len(strains.keys()))

###############################################################################
###############################################################################

total_alignments = 0
perfect_alignment = defaultdict(Counter())
alignments = defaultdict(Counter())
num_ignore_aln = 0
prev_read_obj = Read(None)
# sorted_bamfile_name = sort_by_name(bamfile_name,sorted_prefix=prefix+".nameSorted", cores=args['threads'])

# Read the BAM file
bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
for aln in bamfile.fetch(until_eof=True):
    total_alignments += 1

    read = Read(aln)
    if prev_read_obj.name == '':
       prev_read_obj = read

    if prev_read_obj.name != read.name:
        vote_value = prev_read_obj.tally_votes()
        print(prev_read_obj.name+'\t:\t'+vote_value)
        prev_read_obj = read

    # Also test with PCR duplicates removed as well because it could inflate the read proportion by strain. 
    # Since I'm not looking at coverage, just depth.
    if read.keep_alignment():
        num_perfect_alignments += 1
    else:
        num_ignore_aln += 1

    # For the last read
    vote_value = prev_read_obj.tally_votes()
    print(prev_read_obj.name+'\t:\t'+vote_value)

###############################################################################
###############################################################################

Total_Reads_Aligned = len(read_fractions['Aligned'].keys())
logging.info('Read %d alignments for %d reads', num_aln, Total_Reads_Aligned)

# for frac in read_fractions.keys():
#     logging.info('%s fraction of all alignments : %7.3f%%', frac, len(read_fractions[frac].keys())*100/Total_Reads_Aligned)

strain_primary_abundance = Counter()
primary_strain_weights = Counter()
perfect_reads_considered = Counter()
for read in perfect.keys():
    if not escrow[read]:
    # Spread the 1 vote/read evenly amongst all the strains with primary alignment for this read
    ## Should always be 1 strain
        perfect_reads_considered[read] += 1
        read_vote_value = 1/len(perfect[read].keys())
        for strain in perfect[read].keys():
            strain_primary_abundance[strain]+=read_vote_value
    
for strain in strain_primary_abundance.keys():
    primary_strain_weights[strain] = strain_primary_abundance[strain]/Total_Reads_Aligned

# primary_strain_weights        
strain_escrow_abundance = Counter()
escrow_strain_weights = Counter()
for read in escrow.keys():
    if perfect[read]:
        escrow_reads_considered[read] += 1
        for strain in escrow[read].keys():
            # Spread 1 vote/read based on the primary weight of strains aligned to
            read_vote_value = 1*primary_strain_weights[strain]
            strain_escrow_abundance[strain] += read_vote_value 
    else:
        pass
            
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