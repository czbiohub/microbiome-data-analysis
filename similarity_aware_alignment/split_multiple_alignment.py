import argparse, os, pysam
import numpy as np
from subprocess import call


def split_multiple_alignment(bamfile_name, max_alignments):
    """
    bamfile_name: full path to bamfile name sorted by name (input)
    entry_type: currently can be unique, multiple_all, multiple_within_genome, multiple_across_genome
                description: goes through the bamfile once and generates all 3 bam files, containing reads that
                align uniquely to one genome, reads that align multiple times to the same genome, and
                reads that align multiple times to multiple genomes. Then downstream functions can
                take one output bam or multiple output bam, create a depth file, and tabulate coverage.
    max_alignments: max number of alignments, delete these reads

    revision: 
    2018.11.08: handling both read 1 and read 2; add bam header
    2018.12.08: made into a separate function instead of inside a class.
                used for spliting a bam file into 3 bam files.
    """
    # Open the input bam file
    bamfile = pysam.AlignmentFile(bamfile_name, mode='rb')
    
    # Create 3 file names for unique_alignments, multiple_alignments_unique_genome, 
    # and multiple_alignments_multiple_genomes
    file_root = bamfile_name.rstrip('.bam')
    unique_alignment_filename = file_root + '.unique_alignments.bam'
    multiple_alignments_unique_genome_filename = file_root + '.multiple_alignments_unique_genome.bam'
    multiple_alignments_multiple_genomes_filename = file_root + '.multiple_alignments_multiple_genomes.bam'
    bam_root = file_root
    
    # Open all three files
    # Use the template argument to copy of the bam header from another file
    with pysam.AlignmentFile(unique_alignment_filename, "wb", template=bamfile) as f1, pysam.AlignmentFile(multiple_alignments_unique_genome_filename, "wb", template=bamfile) as f2, pysam.AlignmentFile(multiple_alignments_multiple_genomes_filename, "wb", template=bamfile) as f3:

        # Treating all reads as having the potential to have multiple alignments
        # tempCov = pd.DataFrame(index=coverage.index)
        # tempCov[arguments.sample_name] = 0
        reference_mapped_to = [] # temp storage of all the strain genomes mapped to
        temp_alignments = [] # Used to temporarily store alignments
        current_read_name = ''
        counter = 0

        # Alignments should be sorted by name, but it is still likely that alignments of the 
        # same read pair are interleaved.
        for alignment in bamfile.fetch(until_eof=True): # until_eof is important because it maintains order of alignment within file
            # removed alignment.is_read1, file should still be sorted by name. 2018.12.08 REMOVE NON CONCORDANT ALIGNMENTS
            if alignment.is_paired and alignment.is_proper_pair and ((alignment.is_reverse and not alignment.mate_is_reverse) or (not alignment.is_reverse and alignment.mate_is_reverse)) and len(alignment.cigartuples) == 1:
                # If it's a new read and If there is at least 1 alignment, aka the array is not empty
                if alignment.query_name != current_read_name:
                    # Must have this as a separate if because otherwise reference_mapped_to will never populate
                    if reference_mapped_to:
                        # If the read aligns uniquely
                        # be careful here because the two alignments might be from different read pairings (unlikely though)
                        # 2018.12.08 IGNORE ORPHAN PAIRS HERE. ONLY RELEVANT FOR UNIQUE MAPPED READS
                        if len(reference_mapped_to) == 2:
                            # print('unique alignment')
                            # Output the alignment into the unique alignment file
                            for a in temp_alignments:
                                t = f1.write(a)
                        # else read aligns to multiple locations
                        # If the number of alignments exceed max_alignments, don't process aka ignore read
                        elif len(reference_mapped_to) > 2 and len(reference_mapped_to) < max_alignments:
                            # De-replicate the reference list
                            # For reads that map to multiple references, I'm not removing those that map to off limit regions
                            reference_mapped_to = list(set(reference_mapped_to)) # unique elements
                            # If the read aligns to multiple locations but on the same genome
                            if len(reference_mapped_to) == 1:
                                # print('multiple alignment to single genome')
                                # Output all alignments to multiple but single genome file
                                for a in temp_alignments:
                                    t = f2.write(a)
                            # else the read aligns to multiple genomes
                            else:
                                # print('multiple alignments to multiple genomes')
                                # Output all alignments to multiple genomes file
                                for a in temp_alignments:
                                    t = f3.write(a)
                    # Update the variables to handle the next read
                    reference_mapped_to = []
                    temp_alignments = []
                    current_read_name = alignment.query_name
                # If still processing the same read pair, continue
                ref = bamfile.get_reference_name(alignment.reference_id)
                reference_mapped_to.append(ref.split('_')[0])
                temp_alignments.append(alignment) # append new alignment
            # Counter to maintain progress, for all alignments, not just proper ones
            counter += 1
            if counter % 100000 == 0:
                print('.', end='')
                counter = 0
    bamfile.close()
    print('\n')

        
# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python split_multiple_alignment.py input_bamfile max_alignment_number"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='bamfile', action='store', type=str)
    p.add_argument(dest='maxAlignment', action='store', type=int)

    A = p.parse_args()

    try:
        split_multiple_alignment(A.bamfile, A.maxAlignment)

    except ValueError as e:
        print("ERROR: ValueError:",e)
        print(usage)
    except TypeError as e:
        print("ERROR: TypeError:",e)
        print(usage)
    except IOError as e:
        print("ERROR: IOError %s",e)
        print(usage)

