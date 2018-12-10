import argparse, os, pysam
import numpy as np
import pandas as pd
from subprocess import call
from multiprocessing import Process, Queue
from time import sleep


class tabulate_alignment_fragment:
    """
    tabulate_alignment_fragment.py
    
    reference_list_file: File with two columns used to generate bowtie2 index.
                         Stored with the bowtie2 references and called genome_name_list.csv
    unique_alignment_file: bam file with all unique alignments
    multiple_alignment_file: bam file with all reads that mapped multiple times
    tabulated_alignment_file: output file including all the genomes and one column with fragments
    sample_name: a string that represents the name of the sample
    
    Description: This function is used in accu_align_v1.sh. Only the first read alignments are used.
                 For unique alignments, the number of bps in inferred fragment length is used.
                 For multiple alignments, the number of bps are proportionally assigned.
                 
    Author: Brian Yu
    
    Revision History:
    2018.09.10 Version 1 had the wrong algorithm. This version uses the unique reads.
    2018.12.08 Looking at both paired alignments, ignore non-concordant alignments,
               remove nonperfect CIGAR, remove orphan reads that are unique. 
               Currently does not remove multi-mapped reads that are orphans.
    """

    def  __init__(self, reference_list_file, sample_name):
        """
        reference_list_file: arguments.reference_list_file
        sample_name: arguments.sample_name
        """
        self.sample_name = sample_name
        # Read in genome names, use the corrected name as the row indices
        reference_list = pd.read_csv(reference_list_file, sep=',', header=0, index_col=1)
        reference_list[self.sample_name] = 0 # Added a column with 0's with sample_name as header
        self.coverage = reference_list.drop('file_name', axis=1) # file_name is column header, axis=1 means column

    def process_unique_alignments(self, bamfile_name):
        """
        bamfile_name: arguments.unique_alignment_file
        """
        # Process unique alignments
        print("Processing Uniquely Aligned Reads")
        bamfile = pysam.AlignmentFile(bamfile_name, mode='rb')
        # alignment here is a pysam AlignedSegment data structure
        counter = 0
        for alignment in bamfile.fetch(until_eof=True):
            if alignment.is_paired and alignment.is_read1 and alignment.is_proper_pair and ((alignment.is_reverse and not alignment.mate_is_reverse) or (not alignment.is_reverse and alignment.mate_is_reverse)):
                ref = bamfile.get_reference_name(alignment.reference_id)
                self.coverage.ix[ref.split('_')[0], self.sample_name] += float(abs(alignment.template_length))
                counter += 1
                if counter % 100000 == 0:
                    print('.', end='')
                    counter = 0
        bamfile.close()
        print('.')
        # For debugging purposes
        # self.coverage.to_csv(arguments.tabulated_alignment_file, index=True, header=True)


    def process_multiple_alignments(self, bamfile_name):
        """
        bamfile_name: arguments.multiple_alignment_file
        CIGAR string is handled previously when spliting the bam. It is not handled here.
        """

        # Process multiple alignments
        print("Processing Reads With Multiple Alignments")
        # Create another df with same index
        tempCov = pd.DataFrame(index=self.coverage.index)
        tempCov[self.sample_name] = 0
        reference_mapped_to = [] # temp storage of all the strain genomes mapped to
        inferred_template_size = [] # temp storage of all the inferred template size
        current_read_name = ''
        # Read in bamfile
        bamfile = pysam.AlignmentFile(bamfile_name, mode='rb')
        counter = 0

        for alignment in bamfile.fetch(until_eof=True):
            # if alignment.is_paired and alignment.is_read1 and alignment.is_proper_pair:
            # this part is still using all reads even if half of it falls into an off limit region
            if alignment.is_paired and alignment.is_proper_pair and ((alignment.is_reverse and not alignment.mate_is_reverse) or (not alignment.is_reverse and alignment.mate_is_reverse)):
                # If it's a new read
                if alignment.query_name != current_read_name:
                    # Split the read and add to the original coverage df
                    if reference_mapped_to:
                        reference_mapped_to = list(set(reference_mapped_to)) # unique elements
                        tempSeries = self.coverage.loc[reference_mapped_to, self.sample_name]
                        # If, based on the unique reads, there are more than 0 total bps across ref mapped to
                        if tempSeries.sum() > 0:
                            tempCov.loc[reference_mapped_to, self.sample_name] += tempSeries / tempSeries.sum() * np.median(inferred_template_size)
                        # else: # Do not process reads like this
                            # print(alignment.query_name+' split evenly between '+str(tempSeries.size)+' genomes.')
                            # tempCov.loc[reference_mapped_to, arguments.sample_name] += np.median(inferred_template_size) / tempSeries.size
                    # Update the variables to handle the next read
                    reference_mapped_to = []
                    inferred_template_size = []
                    current_read_name = alignment.query_name
                ref = bamfile.get_reference_name(alignment.reference_id)
                reference_mapped_to.append(ref.split('_')[0])
                inferred_template_size.append(float(abs(alignment.template_length)))
                counter += 1
                if counter % 100000 == 0:
                    print('.', end='')
                    counter = 0
        bamfile.close()
        # For debugging
        # print(coverage)
        # print(tempCov)
        self.coverage.loc[self.coverage.index, self.sample_name] += tempCov.loc[tempCov.index, self.sample_name]
        print('.')


    def output_coverage(self, output_file_name):
        """
        output_file_name: arguments.tabulated_alignment_file
        """
        # Output final tabulated file
        self.coverage.to_csv(output_file_name, index=True, header=True)



class tabulate_basepair_coverage:
    """
    Class Description:  Create coverage files across the genome of all reference sequences. Use 3 bam
                        files - unique; multiple-unique; and multiple-multiple
    Revision History:   2018.10.11 Brian Yu Created. Need to think about how to use multiprocessing
                        2018.10.13 Use df.groupby in compute coverage
                        2018.10.14 Completed first path and completed debugging
                        2018.11.10 Resolved all bugs. First working version of tabulation pipeline
                        2018.12.09 Removed the function to split bam files because it's included 
                                   in a separate function. 
                                   Integrated into the similarity aware alignment method
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
        Description: If depth file is empty, nothing gets produced.
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
            # But the file might be empty, so you need to check if coverage_block is empty
            if bool(previous_contig_name):
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

def DNA_count(ref, sample_name, uu, mm, mu, output):
    f = tabulate_alignment_fragment(ref, sample_name)
    f.process_unique_alignments(uu)
    f.process_multiple_alignments(mm)
    f.process_multiple_alignments(mu)
    f.output_coverage(output)

def DNA_coverage(ref, window_size, uu, mu, mm):
    k = tabulate_basepair_coverage(ref, []) # use default [] argument for ref_list
    print(uu + '\n' + mu + '\n' + mm)
    # unique alignment
    bam_root = '.'.join(uu.split('.')[0:-2])
    k.output_depth_file(bam_root+'.depth_unique.txt', k.create_depth_file(uu, []))
    k.extract_lines_from_depth_file(bam_root+'.depth_unique.txt', 0, window_size, bam_root+'.coverage_unique.csv')
    # multiple to unique
    bam_root = '.'.join(mu.split('.')[0:-2])
    k.output_depth_file(bam_root+'.depth_multiple_unique.txt', k.create_depth_file(mu, []))
    k.extract_lines_from_depth_file(bam_root+'.depth_multiple_unique.txt', 0, window_size, bam_root+'.coverage_unique_multiple.csv')
    # multiple to multiple
    bam_root = '.'.join(mm.split('.')[0:-2])
    k.output_depth_file(bam_root+'.depth_multiple_multiple.txt', k.create_depth_file(mm, []))
    k.extract_lines_from_depth_file(bam_root+'.depth_multiple_multiple.txt', 0, window_size, bam_root+'.coverage_multiple_multiple.csv')



# When running the script from command line, the following lines are executed
if __name__ == "__main__":
    usage = "USAGE: python tabulate_alignment_fragment.py [options] reference_fasta reference_list unique multiple_to_multiple multiple_to_unique output_file"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    p.add_argument('-s', dest='sample_name', action='store', type=str, default='test_sample')
    p.add_argument('-w', dest='window_size', action='store', type=int, default=1000)
    p.add_argument(dest='ref_fasta', action='store', type=str)
    p.add_argument(dest='reference_list_file', action='store', type=str)
    p.add_argument(dest='unique_alignment_file', action='store', type=str)
    p.add_argument(dest='multiple_align_to_multiple', action='store', type=str)
    p.add_argument(dest='multiple_align_to_unique', action='store', type=str)
    p.add_argument(dest='tabulated_alignment_file', action='store', type=str)

    A = p.parse_args()

    try:
        p1 = Process(target=DNA_count, args=(A.reference_list_file, A.sample_name, A.unique_alignment_file, A.multiple_align_to_multiple, A.multiple_align_to_unique, A.tabulated_alignment_file))
        p1.start()
        p1.join()
        p2 = Process(target=DNA_coverage, args=(A.ref_fasta, A.window_size, A.unique_alignment_file, A.multiple_align_to_unique, A.multiple_align_to_multiple))
        p2.start()
        p2.join()
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


