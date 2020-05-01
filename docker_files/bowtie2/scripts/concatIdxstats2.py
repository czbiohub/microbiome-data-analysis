#!/usr/bin/env python3

import logging
import subprocess
import sys
from random import sample

import pandas as pd
import s3fs


def get_idxfiles_list():
    my_bucket = 'czbiohub-microbiome'
    tmp_all_filepaths = subprocess.run("aws s3 ls s3://czbiohub-microbiome/Sonnenburg_Lab/Project_Vital/Bowtie2/Alignments/ --recursive | awk '{print $4}' | grep 'idxstats.txt'",
                                    shell = True,
                                    stdout=subprocess.PIPE)
    tmp_all_filenames = str(tmp_all_filepaths.stdout.decode('utf-8')).split('\n')
    del tmp_all_filenames[-1]
    return [f's3://{my_bucket}/{f}' for f in tmp_all_filenames]

def extract_sample_name(idxFile):
    return str(idxFile).split('/')[7]

def read_idxstats(idxFile):
    return pd.read_table(idxFile,
                         names = ['seq_names','seq_len','mapped_reads','unmapped_reads'],
                         header = None,
                         skipfooter = 1,
                         index_col = 'seq_names')

def get_mapped_reads(idxFile, sample_name):
    return read_idxstats(idxFile)\
        .drop(columns = ['seq_len','unmapped_reads'])\
            .rename(columns={"mapped_reads": sample_name})

def get_seq_len_from_idxstats(idxFile):
    return read_idxstats(idxFile)\
        .drop(columns = ['mapped_reads','unmapped_reads'])

logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s\t[%(levelname)s]:\t%(message)s')

logging.info('Started')

all_filenames = get_idxfiles_list()
total_files = len(all_filenames)
logging.info(f'# of files found: {total_files}')

logging.info(f'Extracting db sequence lengths ...')
get_seq_len_from_idxstats(all_filenames[0]).to_csv("gene_cluster.seqLens.csv.gz")

# Process the first dataframe
logging.info(f'Generating initial dataframe ...')
i = 1
first_sample = extract_sample_name(all_filenames[0])
logging.info(f'\t[{i}/{total_files}] Processing: {first_sample}')
combined_csv = get_mapped_reads(all_filenames[0], first_sample)

# combine all files in the list
logging.info(f'Adding columns to the initial dataframe ...')
for idxFile in all_filenames[1:]:
    i += 1
    try:
        sample_name = get_sample_name(idxFile)
    except:
        logging.critical(f'[FATAL] Could not parse sample name with: {idxFile}')
        sys.exit(1)

    logging.info(f'\t[{i}/{total_files}] Processing: {sample_name}')
    combined_csv[sample_name] = get_mapped_reads(idxFile, sample_name)[sample_name]

logging.info(f'Writing a dataframe of {combined_csv.shape} rows and columns.')
# export to csv
combined_csv.to_csv("all_samples_to_gene_cluster.readAlnStats.csv.gz")
logging.info(f'Done!')
