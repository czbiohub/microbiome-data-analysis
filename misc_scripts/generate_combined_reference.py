"""
reference_list: a file containing a list of input and output genome name
                Needs to be in the same input folder as other fasta files
                First column is input file name; second column is output filename (genome name)
input_folder:   The input folder needs to contain only fasta reference files
output_fasta:   The name of the output combined fasta file in the same folder as input_folder

Description:    This function will not combine contigs into one sequence.
                Instead, it changes the names of the contigs to >genome_name_node1
                where genome_name is also the file name (aka desired reference name)
                The function also tabulates genome size and number of contigs.
"""

import argparse
import numpy as np
import pandas as pd

usage = "USAGE: python generate_combined_reference.py reference_list input_folder output_fasta"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument(dest='reference_list_file', action='store', type=str)
p.add_argument(dest='input_folder', action='store', type=str)
p.add_argument(dest='output_fasta_name', action='store', type=str)

arguments = p.parse_args()

# Read in the list of reference names. key is file name, first column is fixed file name
reference_list = pd.read_csv(arguments.input_folder+'/'+arguments.reference_list_file, sep=',', header=0, index_col=0)
column_name = list(reference_list)

# Setup variables to tabulate reference genome statistics
reference_list['number_of_contigs'] = ""
reference_list['size_of_genome'] = ""

# Open the output reference file name
with open(arguments.input_folder+'/'+arguments.output_fasta_name, 'w') as output_file:

    # For each reference
    for ref in reference_list.index:
        
        temp_genome_size = 0
        temp_contig_number = 0
        
        # Read in the reference file as header and seq
        with open(arguments.input_folder+'/'+ref+'.fna', 'r') as f:
            for l in f:
                # If the line is a contig header
                if '>' in l:
                    t = output_file.write('>' + reference_list.ix[ref,column_name[0]] + '_Node_' + str(temp_contig_number) + '\n')
                    temp_contig_number = temp_contig_number + 1
                # else the line is a sequence
                else:
                    t = output_file.write(l)
                    temp_genome_size = temp_genome_size + len(l.rstrip())
                    
        # For debugging purposes
        print(ref + '\t' + str(temp_contig_number) + '\t' + str(temp_genome_size))
        
        reference_list.ix[ref,'number_of_contigs'] = temp_contig_number
        reference_list.ix[ref,'size_of_genome'] = temp_genome_size
        
    # Output into a combined fasta reference file in input_folder
    reference_list.to_csv(arguments.input_folder+'/reference_genome_statistics.csv', index=True, header=True)
