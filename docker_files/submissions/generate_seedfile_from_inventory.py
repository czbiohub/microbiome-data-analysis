#!/usr/bin/env python

# Search through czb-seqbot fastq inventory for specific sample names and
# generate a seedfile for further analysis
## Also requires pyarrow, fastparquet
## conda install -y pyarrow fastparquet

# USAGE: python generate_seedfile_from_inventory.py <sample_names.list> <s3://inventory_s3path> <source_bucket_name> <output_file.csv>


import logging
import pandas as pd
import boto3
import s3fs
import io
import sys


def generate_seedfile(inventory_df, samples_list, source_bucket):
    seedfile_list = list()
    for sample in samples_list:
        try:
            fwd, rev, *extra = sorted(
                # inventory_df[inventory_df["key"][~(df.columns.str.endswith('Name') | df.columns.str.endswith('Code'))]]
                inventory_df[inventory_df["key"].str.contains(sample, na=False)]["key"]
            )
        except ValueError as v:
            logging.error(f"Could not find sample: '{sample}'.")
            raise Exception(v)

        if len(extra) > 0:
            logging.error(f"Error at sample: {sample}")
            logging.error(f"Captured Forward Fastq: {fwd}")
            logging.error(f"Captured Reverse Fastq: {rev}")
            _ = [logging.error(f"Extra Fastqs: {e}") for e in extra]
            raise Exception(
                "Data does not conform to the expected read pair structure. Possible duplication of sample names."
            )

        seedfile_list.append(
            {
                "sampleName": sample,
                "R1": f"s3://{source_bucket}/{fwd}",
                "R2": f"s3://{source_bucket}/{rev}",
            }
        )

    return pd.DataFrame(seedfile_list)


def read_list_file(list_file):
    with open(list_file) as items_file:
        items_list = [i.rstrip("\n") for i in items_file]

    return items_list


def split_s3uri(s3uri):
    s3path_as_list = s3uri.replace("s3://", "").rstrip("/").split("/")
    bucket = s3path_as_list.pop(0)
    prefix = "/".join(s3path_as_list)

    return bucket, prefix


def read_s3_inventory(bucket_name, file_name):
    client = boto3.resource("s3")
    object = client.Object(bucket_name, file_name)

    buffer = io.BytesIO()
    object.download_fileobj(buffer)
    inventory = pd.read_parquet(buffer)

    # filter to only include compressed fastq files
    inventory_df = inventory[
        (
            inventory["key"].str.endswith(".fastq.gz", na=False)
            | inventory["key"].str.endswith(".fq.gz", na=False)
        )
    ]

    return inventory_df


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    # list of sample names as in the sequencing summary sheet
    # no headers; one sample name per line
    samples_list_file = sys.argv[1]
    # s3 location of the inventory file
    s3_inventory_path = sys.argv[2]
    # bucket for which the inventory is being taken
    seq_bucket = sys.argv[3]
    # Output file name
    output_file = sys.argv[4]

    samples_list = read_list_file(samples_list_file)
    bucket_name, file_name = split_s3uri(s3_inventory_path)

    inventory = read_s3_inventory(bucket_name, file_name)
    samples_df = generate_seedfile(inventory, samples_list, source_bucket=seq_bucket)
    samples_df.to_csv(output_file, index=False)
