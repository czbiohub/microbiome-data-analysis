#!/usr/bin/env python3

import concurrent.futures
import logging
import os
import sys
import shutil
import boto3
import pandas as pd
from tqdm import tqdm


def get_file_names(bucket_name, prefix, suffix="txt"):
    """
    Return a list for the file names in an S3 bucket folder.

    :param bucket: Name of the S3 bucket.
    :param prefix: Only fetch keys that start with this prefix (folder name).
    :param suffix: Only fetch keys that end with this suffix (extension).
    """
    s3_client = boto3.client("s3")
    response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
    objs = response["Contents"]

    while response["IsTruncated"]:
        response = s3_client.list_objects_v2(
            Bucket=bucket_name,
            Prefix=prefix,
            ContinuationToken=response["NextContinuationToken"],
        )
        objs.extend(response["Contents"])

    logging.info(f"Sifting through {len(objs)} files ...")

    shortlisted_files = list()
    if suffix == "":
        shortlisted_files = [obj["Key"] for obj in objs]
        total_size_bytes = sum([obj["Size"] for obj in objs])
    else:
        shortlisted_files = [obj["Key"] for obj in objs if obj["Key"].endswith(suffix)]
        total_size_bytes = sum(
            [obj["Size"] for obj in objs if obj["Key"].endswith(suffix)]
        )

    logging.info(
        f"Found {len(shortlisted_files)} files, totalling about {total_size_bytes/1e9:,.3f} Gb."
    )
    return shortlisted_files


def get_mapped_reads(idxFile):
    return pd.read_table(
        idxFile,
        names=["seq_names", "seq_len", "mapped_reads", "unmapped_reads"],
        header=None,
        skipfooter=1,
        index_col="seq_names",
        engine="python",
    )


def read_idxstats(idxFile, output_path):
    sample_name = os.path.basename(idxFile).split("_vs_")[0]
    # logging.info(f"\tProcessing {sample_name} ...")
    tmpfile = os.path.join(output_path, f"{sample_name}.csv")
    # if file exists, return file name, else create file
    if not os.path.isfile(tmpfile):
        (
            read_idxstats(idxFile)
            .drop(columns=["seq_len", "unmapped_reads"])
            .rename(columns={"mapped_reads": sample_name})
            .to_csv(tmpfile)
        )
    return tmpfile


def aggregate_stats(all_file_paths, tmpdirname):
    list_of_csv = list()
    total_files = len(all_file_paths)
    logging.info(f"Extracting relevant data from {total_files} files ...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        # with concurrent.futures.ThreadPoolExecutor(max_workers=None) as executor:
        future = [
            executor.submit(read_idxstats, idxFile, tmpdirname)
            for idxFile in all_file_paths
        ]
        for f in tqdm(
            concurrent.futures.as_completed(future),
            total=total_files,
            ascii=True,
            desc="Extracting",
        ):
            list_of_csv.append(f.result())
            # yield f.result()
        # processed_files = 0
        # for f in concurrent.futures.as_completed(future):
        #     processed_files += 1
        #     logging.info(f"\tProcessed {processed_files}/{total_files}")
        #     # list_of_df.append(f.result())
        #     yield f.result()

    # combined_csv[sample_name] = get_mapped_reads(idxFile, sample_name)[sample_name]
    # export to csv
    return list_of_csv


# def get_seq_len_from_idxstats(idxFile):
#     return read_idxstats(idxFile).drop(columns=["mapped_reads", "unmapped_reads"])


def combine_csv(list_of_csv, prefix):
    total_files = len(list_of_csv)
    logging.info(f"Aggregating data from {total_files} files")
    list_of_dfs = (
        pd.read_csv(csv, index_col="seq_names")
        for csv in tqdm(list_of_csv, ascii=True, total=total_files, desc="Aggregating",)
    )
    pd.concat(list_of_dfs, axis=1, sort=False).to_csv(f"{prefix}.combined.csv")


if __name__ == "__main__":
    logging.basicConfig(
        # filename=logfile,
        # filemode='w+',
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    cores = 10
    bucket = "czbiohub-microbiome"
    prefix = "Sonnenburg_Lab/InfantMicrobiome/Gene_Clusters_Alignment/Alignment"
    suffix = ".idxstats.txt"
    output_prefix = "InfantMicrobiome"

    # Get all S3 paths with the specified extension
    all_file_paths = [
        f"s3://{bucket}/{f}" for f in get_file_names(bucket, prefix, suffix)
    ]

    # Create a temporary directory
    tmpdir = f"tmp__{output_prefix}"
    os.makedirs(tmpdir, exist_ok=True)

    # Download, clean and save the files
    csv_list = aggregate_stats(all_file_paths, tmpdir)

    # Combine all the files into a single csv file
    combine_csv(csv_list, output_prefix)

    # Remove the temporary files and directory
    try:
        shutil.rmtree(tmpdir)
    except OSError as e:
        logging.critical("Error: %s - %s." % (e.filename, e.strerror))

    # Rejoice!
    logging.info("Huzzah!")
