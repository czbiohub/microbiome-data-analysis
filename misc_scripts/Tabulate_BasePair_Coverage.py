import argparse, os, pysam
import numpy as np
import pandas as pd
from subprocess import call
from multiprocessing import Process, Queue
from time import sleep


class Tabulate_BasePair_Coverage:
    """
    Class Description:  
    Revision History:   2018.10.11 Brian Yu Created. Need to think about how to use multiprocessing
                        2018.10.13 Use df.groupby in compute coverage
                        2018.10.14 Completed first path and completed debugging
                        2018.11.10 Resolved all bugs. First working version of tabulation pipeline
    """

    def __init__(self, reference_fasta, reference_names):
        """
        reference_fasta: full path to the fasta file
        reference_names: file with reference names interested or []
        """
        # Process reference contigs
        self.contig_list = {} # a dictionary with contig name and length
        temp_contig_name = None
        temp_contig_size = None
        with open(reference_fasta,'r') as f:
            for l in f:
                if '>' in l:
                    if temp_contig_name and temp_contig_size:
                        self.contig_list[temp_contig_name] = temp_contig_size
                    temp_contig_name = l.rstrip()[1:]
                    temp_contig_size = 0
                else:
                    temp_contig_size += len(l.rstrip())
            # The last contig still needs to be processed
            self.contig_list[temp_contig_name] = temp_contig_size

        # Process reference contigs
        if not reference_names:
            self.ref = list(set([x.split('_')[0] for x in list(self.contig_list.keys())]))
        else:
            with open(reference_names, 'r') as f:
                self.ref = []
                for l in f:
                    self.ref.append(l.rstrip())
        
        # Debugging
        print('Number of references is: ', len(self.ref))
        # print(self.contig_list)


    def extract_bam_entries(self, bamfile_name):
        """
        bamfile_name: full path to bamfile name sorted by name
        entry_type: currently can be unique, multiple_all, multiple_within_genome, multiple_across_genome
        description: goes through the bamfile once and generates all 3 bam files, containing reads that
                     align uniquely to one genome, reads that align multiple times to the same genome, and
                     reads that align multiple times to multiple genomes. Then downstream functions can
                     take one output bam or multiple output bam, create a depth file, and tabulate coverage.
        revision: 
        2018.11.08: handling both read 1 and read 2; add bam header
        """
        # Open the input bam file
        bamfile = pysam.AlignmentFile(bamfile_name, mode='rb')
                
        # Create 3 file names for unique_alignments, multiple_alignments_unique_genome, 
        # and multiple_alignments_multiple_genomes
        file_root = bamfile_name.rstrip('.bam')
        unique_alignment_filename = file_root + '.unique_alignments.bam'
        multiple_alignments_unique_genome_filename = file_root + '.multiple_alignments_unique_genome.bam'
        multiple_alignments_multiple_genomes_filename = file_root + '.multiple_alignments_multiple_genomes.bam'
        self.bam_root = file_root

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
                # removed alignment.is_read1, file should still be sorted by name
                if alignment.is_paired and alignment.is_proper_pair:
                    # If it's a new read and If there is at least 1 alignment, aka the array is not empty
                    if alignment.query_name != current_read_name:
                        # Must have this as a separate if because otherwise reference_mapped_to will never populate
                        if reference_mapped_to:
                            # If the read aligns uniquely
                            # be careful here because the two alignments might be from different read pairings (unlikely though)
                            if len(reference_mapped_to) == 2:
                                # print('unique alignment')
                                # Output the alignment into the unique alignment file
                                for a in temp_alignments:
                                    t = f1.write(a)
                            # else read aligns to multiple locations
                            else:
                                # De-replicate the reference list
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


    def create_depth_file(self, bamfile_name, core_num):
        """
        bamfile_name: name of the bamfile to create depth file with
        core_num: core_num is only used in sorting and is optional. Default is 2
        Description: This function calls samtools.
                     You can access the newly created bamfiles through self.bam_root
                     The output from pysam.depth is returned as a string, no longer saved as a file
        """
        if not core_num:
            core_num = 2
        else:
            core_num = int(core_num)
        sorted_bamfile = bamfile_name.split('.bam')[0] + '.sortedByCoord.bam';
        pysam.sort("--threads",str(core_num),"-m","2G","-o",sorted_bamfile,bamfile_name)
        # The depth file is returned as one single string
        return pysam.depth(sorted_bamfile)

        
    def extract_lines_from_depth_file(self, depth_file_name, sample_number, window_size, output_file_name):
        """
        depth_file_name: output of samtools depth, no header, each column other than the first 2 are samples
                    The depth_file should be sorted by contig and also by position in the contig.
        sample_number: which column to use, starting at 0, could not be an array
        output_file: filename and full path of the output csv file
        """
        # Check that sample_number is not an array
        if type(sample_number) != type(int()):
            raise ValueError('Variable sample_number is not an int.')
        
        # If output_file exist, erase it
        if output_file_name and os.path.isfile(output_file_name):
            os.remove(output_file_name)

        # Open depth file. I'm going line by line here to minimize memory usage
        # aka. I don't want to use pd.read_table to read in the entire file
        # And thus, the choice of pulling out only one column to tabulate is also a design choice
        with open(depth_file_name,'r') as depth_file:
            coverage_block = {'contig_position' : [], 'contig_depth' : []} # define empty dictionary with two fields
            previous_contig_name = None
            contig_name = None
            for l in depth_file:
                contig_name = l.split()[0]
                if contig_name.split('_')[0] in self.ref:
                    # If previous_contig_name exists and it is not the same as the new contig
                    if previous_contig_name and contig_name != previous_contig_name:
                        # Process the last contig coverage
                        self.output_one_contig_coverage(self.compute_coverage_for_one_contig(previous_contig_name, coverage_block, window_size), output_file_name)
                        # Update the variables for the next contig coverage
                        coverage_block = {'contig_position' : [], 'contig_depth' : []}
                    # Update previous_contig_name
                    previous_contig_name = contig_name
                    # split fields in the line (delimited by tab)
                    tmp_line = l.split()
                    # Append position and depth. field 1 = contig name, field 2 = contig position, field 3 = depth
                    coverage_block['contig_position'].append(int(tmp_line[1]))
                    coverage_block['contig_depth'].append(float(tmp_line[sample_number + 2])) # + 2 because first 2 columns are name and position 
            # since the last contig would not be written, write it here
            self.output_one_contig_coverage(self.compute_coverage_for_one_contig(previous_contig_name, coverage_block, window_size), output_file_name)


    def compute_coverage_for_one_contig(self, contig_name, coverage_block, window_size):
        """
        contig_name: a string representing contig_name
        coverage_block: dictionary with two fields, contig_position and contig_depth, both are lists
        window_size: size of the interval in base pairs to sum or average over in terms of coverage
        """ 
        # self.contig_list is a dictionary where the keys are contig names and the values are lengths
        contig_len = self.contig_list[contig_name]
        # lower_bound = range(1, contig_len, window_size)
        # upper_bound = lower_bound[1:] + [contig_len]
        # Convert the dictionary coverage_block to a dataframe with 1 column called depth. Keys are indices
        coverage_df = pd.DataFrame.from_dict(coverage_block, orient='columns', dtype=int)
        # Use dataframe functions to group by genome position window_size, average coverage, 
        grouped_coverage = coverage_df.groupby(coverage_df['contig_position'].floordiv(window_size)).mean()
        # Extract coverage_depth and fill in none covered locations with 0; First two are contig_name and length
        coverage_vector = [contig_name, str(contig_len)]
        for i in range(int(np.ceil(contig_len/window_size))):
            if i in grouped_coverage.index:
                coverage_vector.append(str(grouped_coverage.ix[i, 'contig_depth']))
            else:
                coverage_vector.append('0')
        return coverage_vector


    def output_depth_file(self, depth_file_name, alignment_depth):
        """
        Write the variable alignment_depth to depth_file_name
        """
        with open(depth_file_name, 'w') as f:
            t = f.write(alignment_depth)


    def output_one_contig_coverage(self, line, output_file):
        """
        line: line to write as a list, needs to be a list of strings, will be joined by commas
        output_file: full path to output file, should be a csv file
        """
        if not output_file:
            print(','.join(line))
        else:
            with open(output_file,'a') as d:
                t = d.write(','.join(line) + '\n')


# Defining the function that uses Tabulate_BasePair_Coverage class
def analyze_one_sample(reference_fasta, ref_list, bamfile_name, window_size, sample_id):
    """
    Use the Tabulate_BasePair_Coverage class to extract the coverage profile
    for one sample.
    Currently I'm not checking for empty files.
    """
    alignment_class = Tabulate_BasePair_Coverage(reference_fasta, ref_list)
    alignment_class.extract_bam_entries(bamfile_name) # output 3 bam files in the same folder
    print(alignment_class.bam_root)
    alignment_class.output_depth_file(alignment_class.bam_root+'.depth_unique.txt', alignment_class.create_depth_file(alignment_class.bam_root+'.unique_alignments.bam', []))
    alignment_class.extract_lines_from_depth_file(alignment_class.bam_root+'.depth_unique.txt', sample_id, window_size, alignment_class.bam_root+'.coverage_unique.csv')
    alignment_class.output_depth_file(alignment_class.bam_root+'.depth_multiple_unique.txt', alignment_class.create_depth_file(alignment_class.bam_root+'.multiple_alignments_unique_genome.bam', []))
    alignment_class.extract_lines_from_depth_file(alignment_class.bam_root+'.depth_multiple_unique.txt', sample_id, window_size, alignment_class.bam_root+'.coverage_unique_multiple.csv')
    alignment_class.output_depth_file(alignment_class.bam_root+'.depth_multiple_multiple.txt', alignment_class.create_depth_file(alignment_class.bam_root+'.multiple_alignments_multiple_genomes.bam', []))
    alignment_class.extract_lines_from_depth_file(alignment_class.bam_root+'.depth_multiple_multiple.txt', sample_id, window_size, alignment_class.bam_root+'.coverage_multiple_multiple.csv')

        
# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python Tabulate_BasePair_Coverage.py -g ref_list_file -s sample_id reference_fasta window_size input_bamfile_list"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument(dest='ref_fasta', action='store', type=str)
    p.add_argument(dest='window_size', action='store', type=int)
    p.add_argument(dest='input_file_list', action='store', type=str)
    p.add_argument('-g','--genomes',dest='ref_list_file_name', action='store', type=str, default=[])
    p.add_argument('-s','--sample_id',dest='sample_id', action='store', type=int, default=0)

    A = p.parse_args()

    try:
        print(A)
        coreNum = 16
        with open(A.input_file_list, 'r') as f:
            bamQueue = Queue()
            for l in f:
                bamQueue.put(l.rstrip())
        # Divde all the files into groups corresponding to number of cores.
        coreCnt = 0
        proc = []
        while not bamQueue.empty():
            if coreCnt < coreNum:
                p = Process(target=analyze_one_sample, args=(A.ref_fasta, A.ref_list_file_name, bamQueue.get(), A.window_size, A.sample_id))
                p.start()
                proc.append(p)
                coreCnt += 1
                sleep(10) # put 10 sec in between calling new threads 
            else: # all 16 cores are used so wait until things finish
                for p in proc:
                    p.join()
                proc = []
                coreCnt = 0
        # Finish up the currently running processes
        for p in proc:
            p.join()
        # Print completion statement
        print('Script completed')

    except ValueError as e:
        print("ERROR: ValueError:",e)
        print(usage)
    except TypeError as e:
        print("ERROR: TypeError:",e)
        print(usage)
    except IOError as e:
        print("ERROR: IOError %s",e)
        print(usage)

