#!/usr/bin/env python3

import concurrent.futures
import logging
import os
import sys
import shutil
import boto3
import pandas as pd
from tqdm import tqdm

from sqlalchemy import (
    Table,
    Column,
    String,
    Integer,
    create_engine,
    MetaData,
    UniqueConstraint,
    func,
    ForeignKey,
)


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


def get_mapped_reads(idxFile, sample_idx, tmp_file_path):
    # read the dataframe from S3 and remove rows where mapped reads is 0
    (
        pd.read_table(
            idxFile,
            usecols=[0, 2],
            delimiter="\t",
            names=["gene_name", "count"],
            dtype={"gene_name": "str", "count": "int32"},
            header=None,
            skipfooter=1,
            engine="python",
        )
        .rename_axis("gene_id")
        .query(f"count > 0")
        .assign(sample_id=sample_idx)
        .to_csv(tmp_file_path)
    )
    return


def read_idxstats(idxFile, sample_idx, output_path):
    sample_name = sample_name_from_path(idxFile)
    tmp_file_path = os.path.join(output_path, f"{sample_name}.csv")

    # if file exists, return file name, else create file
    if not os.path.isfile(tmp_file_path):
        get_mapped_reads(idxFile, sample_idx, tmp_file_path)

    return tmp_file_path, {sample_name: sample_idx}


def aggregate_stats(all_file_paths, tmpdirname, cores):
    list_of_csv = list()
    total_files = len(all_file_paths)
    sample_idx_dict = dict()
    logging.info(f"Extracting relevant data from {total_files} files ...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        # with concurrent.futures.ThreadPoolExecutor(max_workers=None) as executor:
        future = [
            executor.submit(read_idxstats, idxFile, sample_idx, tmpdirname)
            for sample_idx, idxFile in enumerate(all_file_paths)
        ]
        for f in tqdm(
            concurrent.futures.as_completed(future),
            total=total_files,
            ascii=True,
            desc="Extracting",
        ):
            csv_file_path, sample_dict = f.result()
            list_of_csv.append(csv_file_path)
            sample_idx_dict.update(sample_dict)

    # export to csv
    return list_of_csv, sample_idx_dict


def combine_csv(list_of_csv, prefix):
    total_files = len(list_of_csv)
    _ = [sample_name_from_path(csv_file) for csv_file in list_of_csv]
    logging.info(f"Aggregating data from {total_files} files")
    list_of_dfs = (
        pd.read_csv(csv, index_col="seq_names")
        for csv in tqdm(list_of_csv, ascii=True, total=total_files, desc="Aggregating",)
    )
    pd.concat(list_of_dfs, axis=1, sort=False).to_csv(f"{prefix}.combined.csv")


def sample_name_from_path(file_path):
    return os.path.basename(file_path).split("_vs_")[0]


def initialize_db(db_path):
    engine = create_engine(db_path, pool_pre_ping=True,)
    metadata = MetaData()
    gene_table = Table(
        "genes",
        metadata,
        Column("gene_id", Integer(), primary_key=True, autoincrement=False,),
        Column("gene_name", String(255), nullable=False, index=True),
        Column("gene_len", Integer(), nullable=False),
    )
    samples_table = Table(
        "samples",
        metadata,
        Column("sample_id", Integer(), primary_key=True, autoincrement=False,),
        Column("sample_name", String(255), nullable=False, index=True),
        Column("raw_reads", Integer(), nullable=False),
        Column("trimmed_reads", Integer(), nullable=False),
        Column("unique_reads", Integer(), nullable=False),
        Column("duplicate_reads", Integer(), nullable=False),
    )
    gene_counts_table = Table(
        "gene_counts",
        metadata,
        Column("gene_id", Integer(), ForeignKey("genes.gene_id")),
        Column("sample_id", Integer(), ForeignKey("samples.sample_id")),
        Column("count", Integer(), nullable=False),
    )

    with engine.connect() as conn:
        metadata.create_all(conn, checkfirst=False)

    return engine, gene_table, samples_table, gene_counts_table


def load_genes_table(engine, table, idxFile):
    table_name = table.fullname
    logging.info(f"Loading values into the '{table_name}' table from {idxFile} ...")
    result = pd.read_table(
        idxFile,
        usecols=[0, 1],
        names=["gene_name", "gene_len"],
        dtype={"gene_name": "str", "gene_len": "int32"},
        header=None,
        skipfooter=1,
        skip_blank_lines=True,
        engine="python",
    ).to_sql(
        table_name,
        con=engine,
        index=True,
        index_label="gene_id",
        if_exists="append",
        chunksize=100000,
        dtype={"gene_name": String(255), "gene_len": Integer()},
    )
    logging.info(f"\tData uploaded. Checking integrity ...")
    with engine.begin() as con:
        rowcount = con.execute(func.count(table.columns.gene_id)).scalar()
        logging.info(f"\t'{table_name}' table contains {rowcount} rows")
    return


def load_samples_table(engine, table, samples_metadata_path, sample_idx_dict):
    table_name = table.fullname
    logging.info(
        f"Loading {len(sample_idx_dict.keys())} sample metadata into the '{table_name}' table from {samples_metadata_path} ..."
    )
    (
        pd.read_csv(
            samples_metadata_path,
            usecols=[
                "sample_name",
                "raw_reads",
                "trimmed_reads",
                "unique_reads",
                "duplicate_reads",
            ],
            dtype={
                "sample_name": "str",
                "raw_reads": "int32",
                "trimmed_reads": "int32",
                "unique_reads": "int32",
                "duplicate_reads": "int32",
            },
            skip_blank_lines=True,
        )
        .assign(sample_id=lambda row: row["sample_name"].map(sample_idx_dict))
        .dropna(axis=0, subset=["sample_id"])
        .to_sql(
            table_name,
            con=engine,
            index=False,
            if_exists="append",
            chunksize=500,
            dtype={
                "sample_id": Integer(),
                "sample_name": String(255),
                "raw_reads": Integer(),
                "trimmed_reads": Integer(),
                "unique_reads": Integer(),
                "duplicate_reads": Integer(),
            },
        )
    )
    logging.info(f"\tData uploaded. Checking integrity ...")
    with engine.begin() as con:
        rowcount = con.execute(func.count(table.columns.sample_id)).scalar()
        logging.info(f"\t'{table_name}' table contains {rowcount} rows")
    return


def load_counts_table(engine, table, csv_file):
    table_name = table.fullname
    sample_name = os.path.splitext(os.path.basename(csv_file))[0]
    logging.info(
        f"Loading values from sample '{sample_name}' into the '{table_name}' table ..."
    )
    pd.read_csv(
        csv_file,
        header=0,
        usecols=["gene_id", "sample_id", "count"],
        dtype={"gene_id": "int32", "sample_id": "int32", "count": "int32"},
    ).to_sql(
        table_name,
        con=engine,
        index=False,
        if_exists="append",
        chunksize=100000,
        dtype={"gene_id": Integer(), "sample_id": Integer(), "count": Integer(),},
    )
    return


if __name__ == "__main__":
    logging.basicConfig(
        # filename=logfile,
        # filemode='w+',
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    cores = 60
    bucket = "czbiohub-microbiome"
    prefix = "Sonnenburg_Lab/InfantMicrobiome/Gene_Clusters_Alignment/Alignment"
    suffix = ".idxstats.txt"
    output_prefix = "InfantMicrobiome"
    samples_metadata = "all_samples_qc_for_db.csv"

    # Get all S3 paths with the specified extension
    all_file_paths = [
        f"s3://{bucket}/{f}" for f in get_file_names(bucket, prefix, suffix)
    ]
    # TEST
    # all_file_paths = [
    #     "s3://czbiohub-microbiome/Sonnenburg_Lab/InfantMicrobiome/Gene_Clusters_Alignment/Alignment/Hadza_MoBio_hadza-I-L_L_13_1278/bowtie2/Hadza_MoBio_hadza-I-L_L_13_1278_vs_db_infantMicrobiome_genes.coord_sorted.idxstats.txt",
    #     "s3://czbiohub-microbiome/Sonnenburg_Lab/InfantMicrobiome/Gene_Clusters_Alignment/Alignment/ERS473026/bowtie2/ERS473026_vs_db_infantMicrobiome_genes.coord_sorted.idxstats.txt",
    #     "s3://czbiohub-microbiome/Sonnenburg_Lab/InfantMicrobiome/Gene_Clusters_Alignment/Alignment/SRS1719428/bowtie2/SRS1719428_vs_db_infantMicrobiome_genes.coord_sorted.idxstats.txt",
    #     "s3://czbiohub-microbiome/Sonnenburg_Lab/InfantMicrobiome/Gene_Clusters_Alignment/Alignment/M8063002/bowtie2/M8063002_vs_db_infantMicrobiome_genes.coord_sorted.idxstats.txt",
    #     "s3://czbiohub-microbiome/Sonnenburg_Lab/InfantMicrobiome/Gene_Clusters_Alignment/Alignment/SRR10643322/bowtie2/SRR10643322_vs_db_infantMicrobiome_genes.coord_sorted.idxstats.txt",
    # ]

    # Create a temporary directory
    tmpdir = f"tmp__{output_prefix}"
    os.makedirs(tmpdir, exist_ok=True)
    db_path = f"sqlite:///{output_prefix}.db"

    # Download, clean and save the files
    csv_list, samples_idx = aggregate_stats(all_file_paths, tmpdir, cores)

    ## Initialize the database
    engine, gene_table, samples_table, gene_counts_table = initialize_db(db_path)
    load_samples_table(engine, samples_table, samples_metadata, samples_idx)
    load_genes_table(engine, gene_table, all_file_paths[0])

    [
        load_counts_table(engine, gene_counts_table, csv_file)
        for csv_file in csv_list
        if csv_file is not None
    ]

    logging.info(f"\tCounts data uploaded. Checking integrity ...")
    with engine.begin() as con:
        rowcount = con.execute(func.count(gene_counts_table.columns.gene_id)).scalar()
        logging.info(f"\t'{gene_counts_table.fullname}' table contains {rowcount} rows")

    # Combine all the files into a single csv file
    # combine_csv(csv_list, output_prefix)

    # # Remove the temporary files and directory
    # try:
    #     shutil.rmtree(tmpdir)
    # except OSError as e:
    #     logging.critical("Error: %s - %s." % (e.filename, e.strerror))

    # Rejoice!
    logging.info("Huzzah!")
