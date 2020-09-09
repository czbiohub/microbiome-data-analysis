#!/usr/bin/env python3
import logging
import sys
import boto3
import pandas as pd

def get_file_names(bucket_name, prefix, suffix="txt"):
    """Get a list of s3paths given certain restrictions on prefix and suffix

    Args:
        bucket_name (str): Name of the S3 bucket.
        prefix (str): Only fetch keys that start with this prefix (folder name).
        suffix (str, optional): Only fetch keys that end with this suffix (extension). Defaults to "txt".

    Returns:
        list: all the file names in an S3 bucket folder.
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

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    bucket = sys.argv[1]
    prefix = sys.argv[2]
    output_file = sys.argv[3]
    
    fwd_suffix = "read1_trimmed.fastq.gz"
    rev_suffix = "read2_trimmed.fastq.gz"


    df = pd.DataFrame()
    df["R1"] = sorted([
        f"s3://{bucket}/{f}" for f in get_file_names(bucket,prefix,fwd_suffix)
    ])

    df["R2"] = sorted([
        f"s3://{bucket}/{f}" for f in get_file_names(bucket,prefix,rev_suffix)
    ])

    df["sampleName"] = df["R1"].apply(lambda x: x.split('/')[-3])

    df.to_csv(output_file, index=False, columns=["sampleName","R1","R2"])

