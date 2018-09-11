"""
seedfile: must be a csv file including all the information for fastq1 fastq2 bamOutput etc.
          the headers must match the required environmental variables in accu_align_v1.sh
output_command: the output file generated with aegea batch commands

Description:    This function generates all the commands needed to submit batch jobs 
                on AWS batch using aegea. Input seedfile required.
"""

import argparse
import pandas as pd

usage = "USAGE: python create_alignment_submission_commands_v1.py seedfile output_command"

# Making default argument list structures
p = argparse.ArgumentParser(usage=usage)
p.add_argument(dest='seedfile', action='store', type=str)
p.add_argument(dest='output_command', action='store', type=str)

arguments = p.parse_args()

# Create base command string
base_string = 'aegea batch submit --ecr-image aligner --storage /mnt=500 --memory 16000 --vcpus 8 --command='
command_string1 = '"cd /mnt; git clone https://github.com/czbiohub/microbiome-data-analysis.git; coreNum=4; memPerCore=2G; maxInsert=3000; maxAlignments=200; genomeReferenceDir=/czbiohub-brianyu/Synthetic_Community/Genome_References/Bowtie2Index_090718; '
command_string2 = 'source /mnt/microbiome-data-analysis/similarity_aware_alignment/accu_align_v1.sh"'

# Read in seedfile, column 1 (ie 2) needs to be sampleName
run_samples = pd.read_csv(arguments.seedfile, sep=',', header=0, index_col='sampleName') 
# print(run_samples)

# Create command file
with open(arguments.output_command, 'w') as command_file:
    for sample in run_samples.index:
        sample_specific_file_names = 'fastq1='+run_samples.loc[sample,'fastq1']+'; fastq2='+run_samples.loc[sample,'fastq2']+'; bamOutput='+run_samples.loc[sample,'bamOutput']+'; relativeAbundanceOutput='+run_samples.loc[sample,'relativeAbundanceOutput']+'; readAccountingOutput='+run_samples.loc[sample,'readAccountingOutput']+'; '
        t = command_file.write('echo "Submitting Sample '+sample+'"\n')
        t = command_file.write(base_string+command_string1+'sampleName='+sample+'; '+sample_specific_file_names+command_string2)
        t = command_file.write('\n\nsleep 60\n\n')

