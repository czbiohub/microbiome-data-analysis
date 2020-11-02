#!/usr/bin/env python

import argparse
import gzip
import os
import re
import sys
import time

import pandas as pd
import requests
from Bio import Entrez


## Functions
def usage():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
        usage=argparse.SUPPRESS,
        description="""Description:
    Given a list of genome names and accession IDs (GCF_# or GCA_#) - one per line, 
    - standardize the genome name
    - attempt to select the latest available assembly URL
    - get NCBI metadata tags;
    Output can directly be used with 'create_ninjamap_db.py' script

    Usage: python get_assembly_metadata_from_ncbi.py  -list scv1_1.s3paths.list -out create_ninjamap_db.seedfile.csv -email your@email.com
    """,
        epilog="""Examples:
    python get_assembly_metadata_from_ncbi.py  -list scv1_1.s3paths.list -out create_ninjamap_db.seedfile.csv -email your@email.com
    """,
    )
    # Required
    p.add_argument(
        "-list",
        dest="db_ref_file",
        action="store",
        type=str,
        required=True,
        help='2 columns ("Name" and "Accession") tab-sep or comma-sep file, one genome per line',
    )
    p.add_argument(
        "-out",
        dest="out",
        action="store",
        type=str,
        required=True,
        help="output database file prefix",
    )
    p.add_argument(
        "-email",
        dest="email",
        action="store",
        type=str,
        required=True,
        help="Valid email to register with NCBI",
    )

    return vars(p.parse_args())


def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    from Bio import Entrez

    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


def get_assembly_records(term, email):
    # provide your own mail here
    Entrez.email = email
    handle = Entrez.esearch(db="assembly", term=term, retmax="200")
    record = Entrez.read(handle)
    return record["IdList"]


def guess_latest_assembly(records):
    latest_records = set()
    uid= ''
    summaries_list = (get_assembly_summary(unique_id) for unique_id in records)
    for summary in summaries_list:
        # If this is empty, means this is the latest assembly
        latestAcc = summary["DocumentSummarySet"]["DocumentSummary"][0][
            "LatestAccession"
        ]
        # points to the latest RefSeq or GenBank ID
        lastMajorAcc = summary["DocumentSummarySet"]["DocumentSummary"][0][
            "LastMajorReleaseAccession"
        ]
        # select the appropriate unique ID
        if "GCF_" in lastMajorAcc:
            uid = summary["DocumentSummarySet"]["DocumentSummary"][0]["RsUid"]
        elif "GCA_" in lastMajorAcc:
            uid = summary["DocumentSummarySet"]["DocumentSummary"][0]["GbUid"]

        if latestAcc == "":
            latest_records.add(uid)
    return latest_records


def select_best_record(records, term):
    """
    'There can only be one truth.'
    If multiple versions exist:
        try to pick the latest.
        If multiple latest versions exist:
            the search term is probably not the best. give up.
    """

    num_records = len(records)
    # Find the best unique ID for the search term
    if num_records == 1:
        return records[0]
    elif num_records == 0:
        print(
            f"Term: '{term}' returned {num_records} records, expected at least 1. Ignoring ..."
        )
        return None
    elif num_records > 5:
        print(
            f"Term: '{term}' returned {num_records} records, expected 1, thus it is not specific enough. Ignoring ..."
        )
        return None
    elif num_records > 1:
        print(
            f"Term: '{term}' returned {num_records} records {records}, expected 1, will attempt to guess the latest ..."
        )
        # get the latest
        latest_records = guess_latest_assembly(records)
        num_latest_records = len(latest_records)
        if num_latest_records == 1:
            unique_id = list(latest_records)[0]
            print(f"{unique_id} seems to point to the latest record.")
            return unique_id
        elif num_latest_records > 1:
            print(
                f"Term: '{term}' returned {num_latest_records} records, expected 1, thus it is not specific enough. Ignoring ..."
            )
            return None
        elif num_latest_records == 0:
            print(
                f"Term: '{term}' returned {num_latest_records} records, expected at least 1, please double check the search term and try again. Ignoring ..."
            )
            return None


def get_assemblies(term, email):
    """Get genbank assembly metadata for a given GCF/GCA ID.
    Args:
        term: search term, GCF/GCA Id
        email: user email address
    """

    unique_id = None
    orgName = None
    label = None
    metadata = None
    link = None
    print(f"Processing: {term}")
    records = get_assembly_records(term, email)
    unique_id = select_best_record(records, term)
    if unique_id is None:
        return pd.Series([unique_id, orgName, label, metadata, link])

    # get summary
    summary = get_assembly_summary(unique_id)
    doc_summ = summary["DocumentSummarySet"]["DocumentSummary"][0]
    tags = get_metadata_tags(doc_summ)
    refseq_exclusion_tags = {
        "excluded-from-refseq",
        "suppressed_refseq",
        "contaminated",
        "anomalous",
    }
    if len(refseq_exclusion_tags.intersection(tags)) > 0 and "latest_genbank" in tags:
        link_preference = "GCA"
    else:
        link_preference = "GCF"
    link, label = get_ftp_link(doc_summ, preference=link_preference)
    if link is None:
        print(
            f"Could not parse any FTP link for Term: '{term}' and record '{unique_id}'. Ignoring ..."
        )
        return pd.Series([unique_id, orgName, label, metadata, link])
    orgName = get_organism_name(doc_summ)
    metadata = ";".join(tags)

    return pd.Series([unique_id, orgName, label, metadata, link])


def get_ftp_link(summary, preference):
    """get FTP link"""
    # get RefSeq ftp link
    if preference == "GCF":
        url = summary["FtpPath_RefSeq"].replace("ftp://", "https://")
        if url == "":
            preference = "GCA"
        else:
            label = os.path.basename(url)
            # get the fasta link - change this to get other formats
            link = os.path.join(url, label + "_genomic.fna.gz")
            # print(link)
            if is_valid_url(link):
                return (link, label)

    # get GenBank ftp link
    if preference == "GCA":
        url = summary["FtpPath_GenBank"].replace("ftp://", "https://")
        if not url == "":
            label = os.path.basename(url)
            # get the fasta link - change this to get other formats
            link = os.path.join(url, label + "_genomic.fna.gz")
            if is_valid_url(link):
                return (link, label)

    return (None, None)


def get_metadata_tags(summary):
    """Get select metadata tags describing the assembly"""
    metadata = set()
    try:
        property_dict = dict(summary["AnomalousList"][0])
    except:
        property_dict = {}

    if "Property" in property_dict:
        contaminated = property_dict["Property"]
    else:
        contaminated = None

    #     RefSeq_category # 'na' == None
    refSeqCategory = summary["RefSeq_category"]
    if (refSeqCategory == "") or (refSeqCategory == "na"):
        refSeqCategory = None

    metadata.add(contaminated)
    metadata.add(refSeqCategory)
    metadata.update(summary["ExclFromRefSeq"])
    metadata.update(summary["PropertyList"])

    # remove 'None'
    return sorted([tag for tag in metadata if tag])


# def get_assembly_info(summary):
#     """Get the assembly ID"""
#     assemblyName = summary["AssemblyName"]


# def get_taxa_info(summary):
#     """Get strain and species taxids"""
#     # Not needed now, maybe later
#     #     Taxid # (Strain TaxID)
#     taxID = summary["Taxid"]
#     speciesName = summary["SpeciesName"]
#     speciesTaxID = summary["SpeciesTaxid"]


def is_valid_url(url):
    if url is None:
        return False

    try:
        time.sleep(1)
        http_url = url.replace("ftp://", "https://")
        response = requests.get(f"{http_url}")
        if response.status_code == 200:
            return True
        else:
            return False
    except:
        return False


def get_organism_name(summary):
    """Standardize formal names and punctuations in NCBI organism names"""
    import re

    organism = ""
    try:
        #         if summary['Biosource']['InfraspeciesList'][0]['Sub_type'] == "strain":
        strain_id = summary["Biosource"]["InfraspeciesList"][0]["Sub_value"].strip()
        if strain_id == "":
            raise Exception("StrainID not found")

        organism = summary["SpeciesName"].strip()
        if strain_id not in organism:
            organism = f"{organism} {strain_id}"
    except:
        print(
            "Biosource does not contain strain name, inferring from Organism name instead."
        )
        # Organism # (Strain Name) but remove the genus name in () in the end
        organism = " ".join(summary["Organism"].strip().split(" ")[:-1])
    # replace all punctuations, special characters and spaces with '-'
    name = re.sub("[^A-Za-z0-9]+", "-", organism)
    # remove '-' if it appears at the end or beginning of the string
    return re.sub("\-+$|^\-+", "", name)


def process_raw_seedfile(seedfile):
    try:
        df = pd.read_csv(seedfile).dropna(axis=0, subset=["Name", "Accession"])
        assert "Name" in df.columns
    except:
        df = pd.read_table(seedfile).dropna(axis=0, subset=["Name", "Accession"])
        assert "Name" in df.columns

    assert "Accession" in df.columns

    return df


if __name__ == "__main__":
    # This file needs to have the following 2 columns: "Name" and "Accession"
    # Usually, this file will be provided by the user (scientist/project owner)
    args = usage()
    db_ref_file = args["db_ref_file"]
    email = args["email"]
    output = args["out"]

    genomes_df = pd.DataFrame()

    # Read input
    genomes_df = process_raw_seedfile(db_ref_file)

    # Get genome Information
    genomes_df[
        ["unique_id", "Strain_Name", "label", "metadata", "FTP_link"]
    ] = genomes_df.apply(
        lambda row: get_assemblies(row["Accession"], email),
        axis=1,
        result_type="expand",
    )
    genomes_df.to_csv(output, index=False)
    print("All done! Huzzah!!")
