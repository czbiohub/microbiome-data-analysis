#!/usr/bin/env python3

# Each read gets 1 vote for candidate strains. 
# The more strains it maps to the lower the value of it's vote.
# Note: we are not calculating coverage, just read assignments.
# MAJOR ASSUMPTION 1: We will always have a complete genome for organisms (exact strain) that we seek in the sample.

import argparse
import gzip
import logging
import os
import pysam
import re
import sys

from collections import defaultdict, Counter
from time import perf_counter as timer

start = timer()
###############################################################################
# Functions for making this more robust
###############################################################################
def human_time(time):
    time = abs(time)
    day = time // (24 * 3600)
    time = time % (24 * 3600)
    hour = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    time_str=format('%02d:%02d:%02d:%02d'%(day,hour,minutes,seconds))
    return time_str

# def bam_is_empty(fn):
#     if os.path.getsize(fn) > 1000000:
#         return False

#     bam = pysam.Samfile(fn, check_sq=False)
#     try:
#         bam.next()
#         return False
#     except StopIteration:
#         return True

# def sort_and_index(file_name, sorted_prefix=None, cores=1):
#     """ Sorts and indexes a bam file by coordinates.
#     """
#     if sorted_prefix is None:
#         sorted_prefix = file_name.replace('.bam', '') + '_sorted'

#     sorted_name = sorted_prefix + '.bam'
#     pysam.sort('-@',cores, file_name, sorted_prefix)
#     pysam.index(sorted_name)

#     return pysam.Samfile(sorted_name, 'rb')

# def sort_by_name(file_name, sorted_name=None, cores=1):
#     """ Sorts a bam file by the read name, for paired-end
#     """
#     if sorted_name is None:
#         sorted_name = file_name.replace('.bam', '') + '_namesorted.bam'

#     pysam.sort('-@',cores,'-n', '-o' , sorted_name, file_name)

#     return pysam.Samfile(sorted_name, 'rb')

###############################################################################
# Side Project
###############################################################################

# def predict_alignment():
#     feature_cols = ["align_len","query_len","aln_cov","quality","perc_id","aln_score","mate_score","mismatches","gap_open","gap_ext","is_dup","is_primary","is_supp"]
#     read the model pickle, predict and write outcome.

###############################################################################
# Setup Input and Script Usage
###############################################################################

usage = """
    USAGE: 
    python ninjaMap.py \
-bam input_bamfile \
-bin tab-delimited file with Col1= contig name and Col2=Bin/Strain name \
-out abundance table output
-log logfile.txt
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
p.add_argument('-fasta', dest='fastafile', action='store', type=str, required = True,
                help='name sorted bam file.')
p.add_argument('-bin', dest='binmap', action='store', type=str,
                help='tab-delimited file with Col1= contig name and Col2=Bin/Strain name')
p.add_argument('-out', dest='output', action='store', type=str,
                help='abundance table output')
p.add_argument('-log', dest='logfile', action='store', type=str,
                help='job log and run summary')
p.add_argument('-stats', dest='statsfile', action='store', type=str,
                help='read stats in tab delimited format')
p.add_argument('-threads', dest='threads', action='store', type=str, default=1,
                help='number of threads available for this job and subprocesses')
p.add_argument('-debug', dest='debug', action='store_true', default=False,    
                help='save intermediate false positives bam file')

args = vars(p.parse_args())

# DEBUG
# bamfile_name = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/Bacteroides-coprophilus-DSM-18228/Bacteroides-coprophilus-DSM-18228.processed.bam"
# abundance_output_file = 'B_coprophilius.ninjaMap.v1.abundance.tsv'
bamfile_name = args['bamfile']
fastafile_name = args['fastafile']
prefix = os.path.basename(bamfile_name).split('.')[0]
vote_file = prefix +'.ninjaMap.votes.tsv.gz'

if args['debug']:
    true_strain = prefix
    tmp_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
    false_positives = pysam.AlignmentFile(prefix + '.ninjaMap.false_positives.bam', "wb", template=tmp_bamfile)

if not args['binmap']:
    binmap_file = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/contig_names_bin_map.txt"
else:
    binmap_file = args['binmap']

if not args['output']:
    abundance_output_file = prefix +'.ninjaMap.abundance.tsv'
else:
    abundance_output_file = args['output']

if not args['statsfile']:
    stats_file = prefix +'.ninjaMap.stats.tsv'
else:
    stats_file = args['statsfile']

if not args['logfile']:
    logfile = prefix +'.ninjaMap.log.txt'
else:
    logfile = args['logfile']

logging.basicConfig(filename=logfile, 
    filemode='w+', 
    level=logging.DEBUG,
    format='%(asctime)s\t[%(levelname)s]:\t%(message)s')

logging.info('Started')
###############################################################################
# Parse the Contig to Bin/Strain name map file.
###############################################################################

logging.info('Processing the Bin Map file: %s ...', binmap_file)

all_strains = defaultdict(float)
bins = defaultdict()
with open(binmap_file, "r") as binmap:
    for line in binmap:
        line=line.rstrip()
        contig_name, strain_name = line.split('\t')
        bins[contig_name] = strain_name
        all_strains[strain_name] += 1

strains_list = list(all_strains.keys())

logging.info('\t%d contigs assigned to %d strains', len(bins.keys()), len(all_strains.keys()))

###############################################################################
# Parse the BAM file
###############################################################################
logging.info('Processing the BAM file: %s ...', bamfile_name)

total_reads = defaultdict(float)
perfect_alignment = defaultdict(Counter)
perfect_alignment_objects = defaultdict( list )

best_aln_length = 0
best_aln_score = 0
num_perfect_alignments = 0
# Read the BAM file
bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
for aln in bamfile.fetch(until_eof=True):
    read_name = aln.query_name
    ref_name = aln.reference_name
    strain = bins[ref_name]
    edit_dist = dict(aln.tags)['NM']
    query_len = aln.query_length
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    align_len = aln.get_overlap(ref_start,ref_end)
    # align_len = aln.query_alignment_length

    # template_len = aln.template_length

    # Categorical Variables
    # is_proper_pair = aln.is_proper_pair
    # is_dup = aln.is_duplicate
    # is_supp = aln.is_supplementary
    # is_primary = not(aln.is_secondary)

    # Continuous Variables
    # aln_cov = align_len/float(query_len)
    # quality = aln.mapping_quality
    # aln_score = dict(aln.tags)['AS']
    # mate_score = dict(aln.tags)['YS']
    # mismatches = dict(aln.tags)['XM']
    # gap_open = dict(aln.tags)['XO']
    # gap_ext = dict(aln.tags)['XG']
    # perc_id = 100*(align_len-edit_dist)/float(align_len)

    # Side Project
    # df = some pandas dataframe with the above values
    # predict_alignment(df)

    # if align_len > best_aln_length:
    #     best_aln_length = align_len

    # if aln_score > best_aln_score:
    #     best_aln_score = aln_score

    # if aln.is_read1:
    #     read_name = read_name + '__1'
    # else:
    #     read_name = read_name + '__2'

    # https://www.biostars.org/p/106126/
    if (edit_dist == 0) and (align_len == query_len):
        perfect_alignment[read_name][strain].append(aln)
        num_perfect_alignments += 1

        if args['debug'] and (not strain == true_strain) and (not perfect_alignment_objects[read_name]):
            perfect_alignment_objects[read_name].append(aln)

    total_reads[read_name] += 1

bamfile.close()

Total_Reads_Aligned = len(total_reads.keys())
Total_Perfect_Alignments = len(perfect_alignment.keys())
logging.info('\tUsed %d reads with perfect alignments, out of %d (%7.3f%%).', 
    Total_Perfect_Alignments,
    Total_Reads_Aligned,
    Total_Perfect_Alignments*100/Total_Reads_Aligned
    )

###############################################################################
# Separate the Primary from the Escrow alignments
# Calculate the Strain abundance distribution based on the primary alignments.
###############################################################################
logging.info('Separating the Primary from the Escrow alignments ...')
strain_primary_abundance = defaultdict(float)
primary_strain_weights = defaultdict(float)
primary_reads_considered = defaultdict(float)
escrow = defaultdict(Counter)
read_vote_contribution = defaultdict(Counter)

for read in perfect_alignment.keys():
    read_vote_value = 0
    if len(perfect_alignment[read].keys()) == 1:
        # Should always be 1 strain
        # True hit. This read gets 1 whole vote to assign to a Strain.
        primary_reads_considered[read] += 1
        read_vote_value = 1
        strain = list(perfect_alignment[read].keys())[0] # basically the only one.
        strain_primary_abundance[strain] += read_vote_value
        read_vote_contribution[read][strain] += read_vote_value

        if args['debug'] and not strain == true_strain:
            aln_obj = perfect_alignment_objects[read][0]

            # Create Bam file
            false_positives.write(aln_obj)
    else:
        # This read was able to align to multiple strains without any issues.
        # Add it to the escrow dict
        for strain in perfect_alignment[read].keys():
            escrow[read][strain] += 1

tmp_bamfile.close()
false_positives.close()

Total_Primary_Reads = len(primary_reads_considered.keys())
logging.info('\tUsed %d reads for primary distribution, out of %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total.', 
    Total_Primary_Reads,
    Total_Perfect_Alignments,
    Total_Primary_Reads*100/Total_Perfect_Alignments,
    Total_Primary_Reads*100/Total_Reads_Aligned,
    )

logging.info('Calculating strain abundance distribution based on Primary alignments ...')
###############################################################################
# Use the strain abundance distribution based on primary alignments to weight
# the Escrow abundance assignments.
###############################################################################

logging.info('Assigning escrow reads based on primary alignment strain abundance ...')

strain_escrow_abundance = defaultdict(float)
escrow_strain_weights = defaultdict(float)
escrow_reads_considered = defaultdict(float)
for read in escrow.keys():
    escrow_reads_considered[read] += 1
    read_vote_value = 0
    normalize_by_primary_wts = 0
    for strain in escrow[read].keys():
        # Normalize by strain weights based on primary alignments
        if strain_primary_abundance[strain]:
            normalize_by_primary_wts += strain_primary_abundance[strain]

    if normalize_by_primary_wts > 0:
        # Spread 1 vote/read based on the primary weight of strains aligned to
        for strain in escrow[read].keys():
            read_vote_value = 1*strain_primary_abundance[strain]/normalize_by_primary_wts
            strain_escrow_abundance[strain] += read_vote_value
            read_vote_contribution[read][strain] += read_vote_value

Total_Escrow_Reads = len(escrow_reads_considered.keys())
logging.info('\tUsed %d reads for escrow distribution, out of %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total.', 
    Total_Escrow_Reads,
    Total_Perfect_Alignments,
    Total_Escrow_Reads*100/Total_Perfect_Alignments,
    Total_Escrow_Reads*100/Total_Reads_Aligned,
    )

###############################################################################
# Check: did any read contribute more than 1 vote?
###############################################################################
votes = gzip.open(vote_file, 'wt')

# Header
line = "Read_Name"
for strain in strains_list:
    line += '\t' + strain

votes.write(line + '\n')

# Votes
line = ''
for read in read_vote_contribution.keys():
    line = read
    for strain in strains_list:
        if read_vote_contribution[read][strain]:
            vote = read_vote_contribution[read][strain]
        else:
            vote=0
        line += '\t' + str(vote)
    votes.write(line + '\n')
    line=''

votes.close()

###############################################################################
# Compute contig coverage
###############################################################################
logging.info('Computing depth and coverage for each strain in the database ...')
cov_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
fasta = pysam.FastaFile(fastafile_name)

bin_meta = defaultdict( lambda: defaultdict(lambda: defaultdict( int )))
for contig_name in bins.keys():
    strain = bins[contig_name]
    for pile in cov_bamfile.pileup(contig = contig_name, stepper = 'all'):
        # Only calculate coverage from reads with perfect alignment.
        # print(":")
        base_cov_contribution = 0
        for read in pile.get_query_names():
            if read_vote_contribution[read]:
                base_cov_contribution += read_vote_contribution[read][strain]

        if base_cov_contribution > 0:
            bin_meta[strain][contig_name]['covered'] += 1
            bin_meta[strain][contig_name]['depth'] += base_cov_contribution
    
    if not bin_meta[strain][contig_name]['covered']:
        bin_meta[strain][contig_name]['covered'] = 0
        bin_meta[strain][contig_name]['depth'] = 0

        
cov_bamfile.close()

coverage = defaultdict( lambda: defaultdict( int ))
for strain in bin_meta.keys():
    genome_size = genome_cover = wt_depth = num_contigs = 0
    
    for contig_name in bin_meta[strain].keys():
        genome_cover += bin_meta[strain][contig_name]['covered']
        genome_size += fasta.get_reference_length(contig_name)
        num_contigs += 1
        if bin_meta[strain][contig_name]['covered'] > 0:
            wt_depth += (bin_meta[strain][contig_name]['depth'] * bin_meta[strain][contig_name]['covered']) / fasta.get_reference_length(contig_name)
        else:
            wt_depth += 0

    coverage[strain]['depth'] = wt_depth/genome_size
    coverage[strain]['percent_cov'] = genome_cover*100/genome_size

    if coverage[strain]['depth'] > 0:
        logging.info('\tCovered %7.3f%% of %s @ a depth of %7.3f from %d contig(s)',
                    coverage[strain]['percent_cov'],
                    strain,
                    coverage[strain]['depth'],
                    num_contigs)
    
###############################################################################
# Aggregate abundance from the Primary and weighted Escrow alignments
# Write to output file
###############################################################################
logging.info('Writing abundance counts to: %s ...', abundance_output_file)

abund = open(abundance_output_file, 'w')
# Header
abund.write('Strain_Name\tAbundance\tWeighted_Depth\tPercent_Covered\n')
for strain in strains_list:
    abundance = (strain_primary_abundance[strain] + strain_escrow_abundance[strain])*100/Total_Reads_Aligned
    # cover = coverage[strain]
    abund.write(strain +'\t'+ 
                str(abundance) +'\t'+ 
                str(coverage[strain]['depth']) +'\t'+ 
                str(coverage[strain]['percent_cov']) +'\n')
abund.close()

stats = open(stats_file, 'w')
# Header
stats.write('File_Name\tReads_Aligned\tReads_wPerfect_Aln\tReads_wPrimary_Votes\tReads_wEscrowed_Votes\n')
stats.write(
    prefix +'\t'+
    str(Total_Reads_Aligned) +'\t'+ 
    str(Total_Perfect_Alignments) +'\t'+
    str(Total_Primary_Reads) +'\t'+
    str(Total_Escrow_Reads) +'\n')
stats.close()

end = timer()
logging.info('Completed in %s (d:h:m:s)', human_time(end - start))