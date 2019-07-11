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
import pandas as pd
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
    time_str = format('%02d:%02d:%02d:%02d'%(day,hour,minutes,seconds))
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
python ninjaMap.py -bin contig_names_bin_map.txt -bam Bacteroides-sp-9-1-42FAA/Bacteroides-sp-9-1-42FAA.processed.sortedByCoord.bam -prefix Bacteroides-sp-9-1-42FAA    
""")
# Required
p.add_argument('-bam', dest='bamfile', action='store', type=str, required = True,
                help='name sorted bam file.')
p.add_argument('-fasta', dest='fastafile', action='store', type=str, required = True,
                help='database fasta file')
p.add_argument('-bin', dest='binmap', action='store', type=str, required = True,
                help='tab-delimited file with Col1= contig name and Col2=Bin/Strain name')
p.add_argument('-outdir', dest='outdir', action='store', type=str, required = True,
                help='output directory')
# Optional
p.add_argument('-prefix', dest='prefix', action='store', type=str,
                help='output prefix')
p.add_argument('-debug', dest='debug', action='store_true', default=False,    
                help='save intermediate false positives bam file')
p.add_argument('-truth', dest='truth', action='store', default=False,    
                help='If using debug, please provide one strain name that you would like to track.')
p.add_argument('-mbq', dest='min_base_qual', action='store', default=20, type=int,    
                help='minimum read base quality to consider for coverage calculations.')

args = vars(p.parse_args())
# DEBUG
# bamfile_name = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/Bacteroides-coprophilus-DSM-18228/Bacteroides-coprophilus-DSM-18228.processed.bam"
# abundance_output_file = 'B_coprophilius.ninjaMap.v1.abundance.tsv'
bamfile_name = args['bamfile']
fastafile_name = args['fastafile']
binmap_file = args['binmap']
output_dir = args['outdir']

os.makedirs(output_dir, exist_ok=True)

if not args['prefix']:
    default_prefix = os.path.basename(bamfile_name).split('.')[0]
    prefix = os.path.join(output_dir, default_prefix)
else:
    prefix = prefix = os.path.join(output_dir, args['prefix'])

abundance_output_file = prefix +'.ninjaMap.abundance.csv'
stats_file = prefix +'.ninjaMap.stats.csv'
vote_file = prefix +'.ninjaMap.votes.csv.gz'
strain_stats_file = prefix +'.ninjaMap.strain_stats.csv'
logfile = prefix +'.ninjaMap.log.txt'

if args['min_base_qual']:
    min_base_qual = args['min_base_qual']
else:
    min_base_qual = 20

if args['debug']:
    true_strain = args['truth']
    tmp_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
    false_positives = pysam.AlignmentFile(prefix + '.ninjaMap.false_positives.bam', "wb", template=tmp_bamfile)
    true_positives = pysam.AlignmentFile(prefix + '.ninjaMap.true_positives.bam', "wb", template=tmp_bamfile)

logging.basicConfig(
    # filename=logfile, 
    # filemode='w+', 
    level=logging.DEBUG,
    format='%(asctime)s\t[%(levelname)s]:\t%(message)s')

logging.info('Started')
###############################################################################
# Classes
###############################################################################
class Strains:
    total_genome_size = 0
    total_strains = 0
    
    def __init__(self, strain_name):
        self.name = strain_name
        self.num_singular_reads = 0
        self.num_escrow_reads = 0
        self.num_contigs = 0
        self.genome_size = 0
        self.cum_primary_votes = 0
        self.cum_escrow_votes = 0
        
        self.total_covered_bases = 0
        self.total_covered_depth = 0
        self.uniquely_covered_bases = 0
        self.uniquely_covered_depth = 0
        self.adj_primary_wt = 0
        
        self.escrow_covered_bases = 0
        self.escrow_covered_depth = 0
        
        self.aln_norm_abundance = 0
        self.genome_norm_abundance = 0
        self.adjusted_votes = 0
        self.escrow_vote_conversion_rate = 0
        self.singular_vote_conversion_rate = 0
        
        self.covered_bases = set()
        self.contigs = defaultdict(int)
        self.singular_bin = defaultdict(int)
        self.escrow_bin = defaultdict(int)

    def __hash__(self):
        return hash(str(self.name))

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)
        
    def add_contig(self, contig_name, contig_length):
        self.contigs[contig_name] = contig_length
        self.genome_size += contig_length
        self.num_contigs += 1
        
    def add_paired_singular_vote(self, read_name, mate_name, read_vote = 1):
#         print(str(self) +'\t:\t' + self.name +'\t:\t' + read.unique_name +'\t:\t'+ str(read_vote))
        # for R1
        self.singular_bin[read_name] += read_vote
        self.cum_primary_votes += read_vote 
        self.num_singular_reads += 1

        # for R2
        self.singular_bin[mate_name] += read_vote
        self.cum_primary_votes += read_vote
        self.num_singular_reads += 1
        
    
    def add_escrow_vote(self, read_name, read_vote):
        if read_name in self.singular_bin.keys():
            return

        self.escrow_bin[read_name] += read_vote
        self.cum_escrow_votes += read_vote
        self.num_escrow_reads += 1
    
    def _add_covered_base(self, contig_name, contig_pos):
        unique_contig_name = contig_name+'_'+str(contig_pos)
        self.covered_bases.add(unique_contig_name)
        
    def calculate_singular_coverage (self, bamfile_name, fasta_file):
        if self.num_singular_reads == 0:
            return
        
        cov_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
        fasta = pysam.FastaFile(fasta_file)
        
        for contig_name in self.contigs.keys():
            for pileupcolumn in cov_bamfile.pileup(contig = contig_name, stepper = 'all', min_base_quality = 20):
                base_cov_contribution = 0
                if pileupcolumn.nsegments > 0:
                    base_cov_contribution = self._calc_base_depth(pileupcolumn, self.singular_bin)

                    if base_cov_contribution > 0:
                        self.uniquely_covered_bases += 1
                        self._add_covered_base(contig_name, pileupcolumn.reference_pos)
                        self.uniquely_covered_depth += base_cov_contribution
        cov_bamfile.close()
        fasta.close()
        
        return
    
    def calculate_escrow_coverage (self,bamfile_name, fasta_file):
        if self.num_escrow_reads == 0:
            return
        
        cov_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
        fasta = pysam.FastaFile(fasta_file)
        for contig_name in self.contigs.keys():
#             for pileupcolumn in cov_bamfile.pileup(contig = contig_name, stepper = 'samtools', fastafile = fasta, min_base_quality = min_base_qual):
            for pileupcolumn in cov_bamfile.pileup(contig = contig_name, stepper = 'all', min_base_quality = 20):
                # For this base/column on the contig
                # print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                base_cov_contribution = 0
                if pileupcolumn.nsegments > 0:
                    base_cov_contribution = self._calc_base_depth(pileupcolumn, self.escrow_bin)

                    if base_cov_contribution > 0:
                        self.escrow_covered_bases += 1
                        self._add_covered_base(contig_name, pileupcolumn.reference_pos)
                        self.escrow_covered_depth += base_cov_contribution
                    
        cov_bamfile.close()
        return

    def _calc_base_depth(self, pileupcolumn, bin_dict):
        wt_base_depth = 0
        for pileupread in pileupcolumn.pileups:
            # For all reads aligned to this base/column on the contig
            # print ('\tbase in read %s = %s' %
            #       (pileupread.alignment.query_name,
            #        pileupread.alignment.query_sequence[pileupread.query_position]))
            read_unique_name = Reads.get_unique_read_name(pileupread.alignment)
            if read_unique_name in bin_dict.keys():
                # Only calculate coverage from reads with perfect alignment(singular or escrow).
                wt_base_depth += bin_dict[read_unique_name]
                self.total_covered_depth += 1

        return wt_base_depth
    
    def _normalize_votes(self):
        # Read to Vote conversion ratio
        self.total_covered_bases = len(self.covered_bases)
        
        if self.num_singular_reads > 0:
            self.singular_vote_conversion_rate = self.cum_primary_votes/self.num_singular_reads
        else:
            self.singular_vote_conversion_rate = 0
            
        if self.num_escrow_reads > 0:
            self.escrow_vote_conversion_rate = self.cum_escrow_votes/self.num_escrow_reads
        else:
            self.escrow_vote_conversion_rate = 0

        if (self.num_singular_reads > 0) or (self.num_escrow_reads > 0):
            self.adjusted_votes = ((self.cum_primary_votes * self.num_singular_reads) + (self.cum_escrow_votes * self.num_escrow_reads))/(self.num_singular_reads + self.num_escrow_reads)
        else:
            self.adjusted_votes = 0
            
        self.aln_norm_abundance = self.adjusted_votes / Reads.reads_w_perfect_alignments
        self.genome_norm_abundance = (self.aln_norm_abundance * self.total_covered_bases) / self.genome_size

            
    
    def compile_general_stats(self):
        '''
        for each object of this class, return a pandas data frame with 
        strains as rows and number of singular and escrow votes
        '''
        if (self.uniquely_covered_bases == 0) and (self.escrow_covered_bases == 0):
            return None
        
        singular_depth = 0
        escrow_depth = 0
        self._normalize_votes()
        
        if self.uniquely_covered_bases > 0:
            singular_depth = (self.uniquely_covered_depth/self.uniquely_covered_bases) * self.genome_size
            
        if self.escrow_covered_bases > 0:
            escrow_depth = (self.escrow_covered_depth/self.escrow_covered_bases) * self.genome_size
        
        read_fraction = (self.cum_escrow_votes+self.cum_primary_votes)*100/Reads.total_reads_aligned
        # Dataframe
        return pd.DataFrame(
            index = [self.name],
            data  = {
                'Genome_Size' : self.genome_size,
                'Total_Bases_Covered' : self.total_covered_bases,
                'Coverage_Depth' : self.total_covered_depth/self.genome_size,
                'Read_Fraction' : read_fraction,
                'Lib_Norm_Abundance' : self.aln_norm_abundance,
                'Genome_Size_Lib_Norm_Abundance' : self.genome_norm_abundance,
                'Adjusted_Votes' : self.adjusted_votes,
                'Singular_Strain_Weight' : self.adj_primary_wt,
                'Total_Singular_Reads' : self.num_singular_reads,
                'Total_Singular_Votes' : self.cum_primary_votes,
                'Singular_Read_Vote_Ratio' : self.singular_vote_conversion_rate,
                'Singular_Fraction_of_all_aligned_Reads' : self.num_singular_reads/Reads.total_reads_aligned,
                'Singular_Coverage' : self.uniquely_covered_bases,
                'Singular_Depth' : singular_depth,
                'Total_Escrow_Reads' : self.num_escrow_reads,
                'Total_Escrow_Votes' : self.cum_escrow_votes,
                'Escrow_Read_Vote_Ratio' : self.escrow_vote_conversion_rate,
                'Fraction_of_all_Escrow_Reads' : self.num_escrow_reads/Reads.total_escrow_reads,
                'Escrowed_Cov' : self.escrow_covered_bases,
                'Escrowed_Depth' : escrow_depth
            }
        )
    
    def compile_by_abundance(self):
        '''
        return a pandas data frame with 1 row x 4 columns. 
                    relative_abundance, percent_coverage, wt_coverage_depth
        strain_name
        '''
        read_fraction = (self.cum_escrow_votes+self.cum_primary_votes)*100/Reads.total_reads_aligned

        if read_fraction == 0:
            return None

        percent_coverage = self.total_covered_bases * 100 / self.genome_size
        coverage_depth = self.total_covered_depth/self.genome_size
        
        # per_base_singular_depth = self.uniquely_covered_depth / self.uniquely_covered_bases
        # per_base_escrow_depth = self.escrow_covered_depth / self.escrow_covered_bases

        wt_coverage_depth = 0
        if (self.uniquely_covered_bases + self.escrow_covered_bases) > 0:
            wt_coverage_depth = (((self.uniquely_covered_depth * self.uniquely_covered_bases) + (self.escrow_covered_depth * self.escrow_covered_bases))/(self.uniquely_covered_bases + self.escrow_covered_bases))/self.genome_size
            # wt_coverage_depth = (per_base_singular_depth + per_base_escrow_depth) * self.total_covered_bases / self.genome_size
        # 'Lib_Norm_Abundance' : 100*self.aln_norm_abundance,
        # 'Genome_Size_Lib_Norm_Abundance' : 100*self.genome_norm_abundance,

        return pd.DataFrame(
                    index = [self.name],
                    data  = {
                        'Read_Fraction' : read_fraction,
                        'Percent_Coverage' : percent_coverage,
                        'Coverage_Depth' : coverage_depth,
                        'Weighted_Coverage_Depth' : wt_coverage_depth
                        }
                    )
    
    def mike_drop(self):
        return self.uniquely_covered_bases / self.genome_size

class Reads:
    total_reads_aligned = 0
    reads_w_perfect_alignments = 0
    total_singular_reads = 0
    total_escrow_reads_kept = 0
    total_escrow_reads_discarded = 0
    total_escrow_reads = 0
    total_singular_reads_in_pairs = 0

    def __init__(self, name):
        self.name = name
        self.unique_name = name
        self.has_voted = False
        self.cum_vote = 0

        self.in_singular_bin = False
        
        self.mates_unique_name = ''
        self.mate_has_perfect_match = False

        self.mapped_strains = defaultdict()

    def __hash__(self):
        return hash(str(self.unique_name))

    def __eq__(self, other):
        return self.unique_name == other.unique_name

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

    def add_exact_match(self, strain):
        self.mapped_strains[strain] = strain.name
    
    def put_pair_in_singular_bin(self, mate):
        self.in_singular_bin = True
        mate.in_singular_bin = True

    def add_vote(self, vote_value):
        self.cum_vote += vote_value
        self.has_voted = True

    def is_fraud(self):
        '''
        for each object of this class, return True if cumulative votes > 1
        '''
        return (round(self.cum_vote, 5) > 1)

    def get_voting_details(self):
        '''
        returns a list of 3 element lists, each containing: strain_name, singular vote value and escrow vote value
        '''
        vote_list = list()
        for strain in self.mapped_strains.keys():
            strain_name = ''
            escrow_votes = 0
            singular_votes = 0
            cumulative_vote = 0
            
            if self.cum_vote is not None:
                cumulative_vote = self.cum_vote

            if strain.name is not None:
                strain_name = strain.name

            if self.unique_name in strain.escrow_bin.keys():
                escrow_votes = strain.escrow_bin[self.unique_name]

            if self.unique_name in strain.singular_bin.keys():
                singular_votes = strain.singular_bin[self.unique_name]

            vote_list.append([strain_name, singular_votes, escrow_votes, cumulative_vote])
        return vote_list
    
    @staticmethod
    def choose_primary_candidate(read, mate):
        if read.mate_has_perfect_match or mate.mate_has_perfect_match :
            read_vote_value = 0
            if read.in_singular_bin and mate.in_singular_bin:
                # they both match a single strain
                read_strain = list(read.mapped_strains.keys())[0] # basically the only one.
                mate_strain = list(mate.mapped_strains.keys())[0] # basically the only one.
                if read_strain.name == mate_strain.name:
                    # the strains are the same
                    return read_strain.name
            elif read.in_singular_bin and not mate.in_singular_bin:
                # R1 matches a single strain, but R2 matches multiple
                read_strain = list(read.mapped_strains.keys())[0] # basically the only one.
                if read_strain.name in mate.mapped_strains.values():
                    # If there is an overlap between strain matches, R1 recruits R2 for it's strain match
                    return read_strain.name
            elif mate.in_singular_bin and not read.in_singular_bin:
                # R2 matches a single strain, but R1 matches multiple
                mate_strain = list(mate.mapped_strains.keys())[0] # basically the only one.
                if mate_strain.name in read.mapped_strains.values():
                    # If there is an overlap between strain matches, R2 recruits R1 for it's strain match
                    return mate_strain.name
        return None

    @staticmethod
    def is_perfect_alignment(aln):
        edit_dist = dict(aln.tags)['NM']
        query_len = aln.query_length
        ref_start = aln.reference_start
        ref_end = aln.reference_end
        
        # https://www.biostars.org/p/106126/
        return ((edit_dist == 0) and (query_len == aln.get_overlap(ref_start, ref_end)))
    
    @staticmethod
    def parse_read_name(aln):
        '''
        Accept: AlignmentFile object from PySam
        if read name has a '/', this is the old format. 
        strip the content after the '/', return remaining
        else, return it as is.
        '''
        try:
            key, value = aln.query_name.split("/")
        except ValueError:
            return str(aln.query_name)
        else:
            return str(key)
        
    @staticmethod
    def get_unique_read_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'fwd'
        else:
            orientation =  'rev'
            
        return Reads.parse_read_name(aln) +'__'+ orientation
    
    @staticmethod
    def get_unique_mate_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'rev'
        else:
            orientation =  'fwd'
            
        return Reads.parse_read_name(aln) +'__'+ orientation

###############################################################################
# Parse the Contig to Bin/Strain name map file.
###############################################################################

logging.info('Processing the Bin Map file: %s ...', binmap_file)

all_strains = defaultdict(list)
with open(binmap_file, "r") as binmap:
    for line in binmap:
        line=line.rstrip()
        contig_name, strain_name = line.split('\t')
        all_strains[strain_name].append(contig_name)

all_strain_obj = defaultdict(int)
strains_list = list(all_strains.keys())
fasta = pysam.FastaFile(fastafile_name)
bins = defaultdict()
for strain_name in strains_list:
    strain = Strains(strain_name)
    all_strain_obj[strain_name] = strain
    for contig_name in all_strains[strain_name]:
        strain.add_contig(contig_name, fasta.get_reference_length(contig_name))
        bins[contig_name] = strain_name
        Strains.total_genome_size += fasta.get_reference_length(contig_name)

logging.info('\t%d contigs assigned to %d strains', len(bins.keys()), len(all_strains.keys()))
fasta.close()

###############################################################################
# Parse the BAM file
###############################################################################
logging.info('Processing the BAM file: %s ...', bamfile_name)

total_reads = set()
if args['debug']:
    perfect_alignment = defaultdict(lambda: defaultdict(list))
else:
    perfect_alignment = defaultdict(lambda: defaultdict(int))
mate_map = defaultdict()

# Read the BAM file
bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
for aln in bamfile.fetch(until_eof=True):
    read_name = Reads.get_unique_read_name(aln)
    mates_name = Reads.get_unique_mate_name(aln)
    mate_map[read_name] = mates_name

    if Reads.is_perfect_alignment(aln):
        strain_name = bins[aln.reference_name] # [contig_name]

        if args['debug']:
            perfect_alignment[read_name][strain_name].append(aln)
        else:
            perfect_alignment[read_name][strain_name] += 1
        
    total_reads.add(read_name)

bamfile.close()

Reads.total_reads_aligned = len(total_reads)
Reads.reads_w_perfect_alignments = len(perfect_alignment.keys())
logging.info('\tUsed %d reads with perfect alignments, out of %d (%7.3f%%).', 
    Reads.reads_w_perfect_alignments,
    Reads.total_reads_aligned,
    Reads.reads_w_perfect_alignments*100/Reads.total_reads_aligned
    )

###############################################################################
# Separate the Primary from the Escrow alignments
# Calculate the Strain abundance distribution based on the primary alignments.
###############################################################################
read_objects = defaultdict()
logging.info('Separating the Primary from the Escrow alignments ...')
for read_name in perfect_alignment.keys():
    read = Reads(read_name)
    read_objects[read_name] = read
    
    mate_name = mate_map[read_name]
    read.mates_unique_name = mate_name
    if mate_name in perfect_alignment.keys():
        # This means the mate had a perfect match too.
        read.mate_has_perfect_match = True

    read.num_strains = len(perfect_alignment[read_name].keys())

    for strain_name in perfect_alignment[read_name].keys():
        read.add_exact_match(all_strain_obj[strain_name])
        if args['debug']:
            if read.num_strains == 1 and strain_name != true_strain:
                for aln in perfect_alignment[read_name][strain_name]:
                    # write to Bam file
                    false_positives.write(aln)
            elif read.num_strains == 1 and strain_name == true_strain:
                for aln in perfect_alignment[read_name][strain_name]:
                    # write to Bam file
                    true_positives.write(aln)

    if read.num_strains == 1:
        # Singular
        # True hit. This read gets 1 whole vote to assign to a Strain.
        read.in_singular_bin = True
        Reads.total_singular_reads += 1

if args['debug']:
    tmp_bamfile.close()
    false_positives.close()
    true_positives.close()

del perfect_alignment
del mate_map

logging.info('\tUsed %d reads for primary distribution, out of %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total.', 
    Reads.total_singular_reads,
    Reads.reads_w_perfect_alignments,
    Reads.total_singular_reads*100/Reads.reads_w_perfect_alignments,
    Reads.total_singular_reads*100/Reads.total_reads_aligned
    )

if len(read_objects.keys()) != Reads.reads_w_perfect_alignments:
    logging.critical('Read %d reads with perfect alignments, but created %d read objects. There is something fishy going on here...',
    Reads.reads_w_perfect_alignments,
    len(read_objects.keys())
    )    
    sys.exit()

votes = gzip.open(vote_file, 'wt')
# Header
votes.write("Read_Name,Strain_Name,cSingular_Vote,cEscrow_Vote,cumulative_votes,was_discarded,is_primary,is_singular\n")

singular_fraud_alert = False
escrow_read_objects = defaultdict()
for name, read in read_objects.items():
    if read.has_voted:
        continue
    
    # if read.mates_unique_name in read_objects.keys():
    if read.mate_has_perfect_match:
        # Means both reads in a mate are prefect alignments
        mate = read_objects[read.mates_unique_name]
        strain_name = Reads.choose_primary_candidate(read, mate)

        if strain_name is not None:
            strain = all_strain_obj[strain_name]
            strain.add_paired_singular_vote(read.unique_name, mate.unique_name, 1)
            read.add_vote(1)
            mate.add_vote(1)
            Reads.total_singular_reads_in_pairs += 2
            votes.write(read.name + ',' + strain.name + ',' + str(strain.singular_bin[read.unique_name]) + ',' + 
                str(strain.escrow_bin[read.unique_name]) + ',' + str(read.cum_vote) + ',False,'+ 
                str(read.has_voted)+','+str(read.in_singular_bin)+'\n')
            
            votes.write(mate.name + ',' + strain.name + ',' + str(strain.singular_bin[mate.unique_name]) + ',' + 
                str(strain.escrow_bin[mate.unique_name]) + ',' + str(mate.cum_vote) + ',False,'+ 
                str(mate.has_voted)+','+str(mate.in_singular_bin)+'\n')
        else:
            escrow_read_objects[name] = read
    else:
        escrow_read_objects[name] = read
    
    if read.is_fraud():
        singular_fraud_alert = True
        num_singular_fraud_reads += 1
    
logging.info('\t%d reads in pairs will be used for singular alignment strain abundance, and %d for escrow alignment strain abundance from %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total.', 
    Reads.total_singular_reads_in_pairs,
    len(escrow_read_objects.keys()),
    Reads.total_singular_reads,
    Reads.total_singular_reads*100/Reads.reads_w_perfect_alignments,
    Reads.total_singular_reads*100/Reads.total_reads_aligned
    )

if len(read_objects.keys()) != (len(escrow_read_objects.keys()) + Reads.total_singular_reads_in_pairs):
    logging.critical('Read %d read objects with perfect alignments, but created %d read objects as singular and escrow total (%d + %d). There is something fishy going on here...',
    len(read_objects.keys()),
    (len(escrow_read_objects.keys()) + Reads.total_singular_reads_in_pairs),
    Reads.total_singular_reads_in_pairs,
    len(escrow_read_objects.keys())
    )    
    sys.exit()

del read_objects
###############################################################################
# Calculate unique number of bases covered for each genome
###############################################################################
logging.info('Computing depth and coverage for each strain in the database based on singular alignments ...')
for name, strain in all_strain_obj.items():
#     print(str(strain) +'\t'+strain.name+'\t'+ str(len(strain.singular_bin.keys())))
    strain.calculate_singular_coverage(bamfile_name, fastafile_name)

###############################################################################
# Use the strain abundance distribution based on Singular alignments to weight
# the Escrow abundance assignments.
###############################################################################
logging.info('Assigning escrow reads based on singular alignment strain abundance ...')
logging.info('\t and making sure there is no voter fraud ...')

escrow_fraud_alert = False
num_fraud_reads = 0

for name, read in escrow_read_objects.items():
    read_vote_value = 0
    total_primary_wts = 0
    to_discard = False

    if not read.has_voted:
        # Allow voting ONLY if the read it hasn't voted already
        # discard = read.calculate_escrow_abundance()
        for strain in read.mapped_strains.keys():
            # Normalize by strain weights based on singular alignments
            if strain.uniquely_covered_bases > 0:
                # The 'Mike adjustment' or The 'Mike drop'
                # strain.adj_primary_wt = strain.uniquely_covered_depth / strain.uniquely_covered_bases
                strain.adj_primary_wt = strain.mike_drop()
                total_primary_wts += strain.adj_primary_wt

        if total_primary_wts > 0:
            Reads.total_escrow_reads_kept += 1
            # Spread 1 vote/read based on the singular weight of strains aligned to
            for strain in read.mapped_strains.keys():
                read_vote_value = 1 * strain.adj_primary_wt / total_primary_wts
                strain.add_escrow_vote(read.unique_name, read_vote_value)
                read.add_vote(read_vote_value)
        else:
            Reads.total_escrow_reads_discarded += 1
            to_discard = True

        # return to_discard
    voting_details = read.get_voting_details()
    for voting_detail_sublist in voting_details:
        strain_name, singular_vote_count, escrow_vote_count, cumulative_vote = voting_detail_sublist
        votes.write(read.name + ',' + strain_name + ',' + str(singular_vote_count) + ',' + str(escrow_vote_count) + ',' + str(cumulative_vote) + ',' + str(to_discard) +','+ str(read.has_voted)+','+str(read.in_singular_bin)+'\n')
        if read.is_fraud():
            escrow_fraud_alert = True
            num_fraud_reads += 1

votes.close()

del escrow_read_objects

Reads.total_escrow_reads = Reads.total_escrow_reads_kept + Reads.total_escrow_reads_discarded

logging.info('\tUsed %d reads for escrow distribution, out of %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total.',
            Reads.total_escrow_reads,
            Reads.reads_w_perfect_alignments,
            Reads.total_escrow_reads*100/Reads.reads_w_perfect_alignments,
            Reads.total_escrow_reads*100/Reads.total_reads_aligned
    )
if Reads.total_escrow_reads_discarded > 0:
    logging.info('\t%d out of %d escrow reads had to be discarded because their singular weight was too low, but %d (%7.3f%%) were still retained',
                Reads.total_escrow_reads_discarded,
                Reads.total_escrow_reads,
                Reads.total_escrow_reads_kept,
                Reads.total_escrow_reads_kept*100/Reads.reads_w_perfect_alignments)
    
if singular_fraud_alert or escrow_fraud_alert:
    logging.critical('[FATAL] There were signs of voter fraud. See the votes file (%s) for the complete picture.', vote_file)
    logging.critical("[FATAL] Voter fraud committed by " +str(num_fraud_reads)+ " out of " +str(Reads.total_escrow_reads)+ " escrow reads")
    sys.exit("Voter Fraud Detected!")

###############################################################################
# Sanity Check
###############################################################################
# for name, strain in all_strain_obj.items():
#     singular_vote_ratio = escrow_vote_ratio = 0
#     if strain.num_singular_reads > 0:
#         singular_vote_ratio = strain.cum_primary_votes/strain.num_singular_reads
#     if strain.num_escrow_reads > 0:
#         escrow_vote_ratio = strain.cum_escrow_votes/strain.num_escrow_reads
                             
#     cumulative_vote = (strain.cum_escrow_votes + strain.cum_primary_votes)/Reads.total_reads_aligned
#     if cumulative_vote > 10e-2:
# #     if strain.num_singular_reads > 0:
#         logfile.info(name + '\t' + str(strain.num_singular_reads) + '\t' + str(strain.num_escrow_reads) + '\t' +  str(singular_vote_ratio) + '\t' + str(escrow_vote_ratio) + '\t' + str(cumulative_vote))

###############################################################################
# Calculate unique number of bases covered for each genome
# Create strain stats file
# Create relative abundance output file
###############################################################################

logging.info('Computing depth and coverage for each strain in the database based on escrow alignments ...')

abundance_df = pd.DataFrame()
stats_df = pd.DataFrame()

for name, strain in all_strain_obj.items():
    strain.calculate_escrow_coverage(bamfile_name, fastafile_name)
    
    strain_stats_df = strain.compile_general_stats()
    if strain_stats_df is not None:
        stats_df = pd.DataFrame.add(stats_df, strain_stats_df, fill_value = 0)

    strain_abundance_df = strain.compile_by_abundance()
    if strain_abundance_df is not None:
        abundance_df = pd.DataFrame.add(abundance_df, strain_abundance_df, fill_value = 0)

logging.info('Writing strain stats file ...')
stats_df.to_csv(strain_stats_file, index_label='Strain_Name')

logging.info('Creating the relative abundance file ...')
abundance_df.to_csv(abundance_output_file, index_label='Strain_Name')

###############################################################################
# Stats
###############################################################################
stats = open(stats_file, 'w')
# Header
stats.write('File_Name,Reads_Aligned,Reads_wPerfect_Aln,Reads_wSingular_Votes,Reads_wEscrowed_Votes,Discarded_Reads_w_Perfect_Aln\n')
stats.write(
    prefix +','+
    str(Reads.total_reads_aligned) +','+ 
    str(Reads.reads_w_perfect_alignments) +','+
    str(Reads.total_singular_reads_in_pairs) +','+
    str(Reads.total_escrow_reads_kept) +','+
    str(Reads.total_escrow_reads_discarded) +'\n')
stats.close()

###############################################################################
# The End
###############################################################################
end = timer()
logging.info('Completed in %s (d:h:m:s)', human_time(end - start))
