#!/usr/bin/env python3
"""
General purpose submission script for pipelines that require you to export S3 paths for:
    Output directory    --> 'S3OUTPUTPATH', 
    Input fwd_fastq     --> 'R1', and 
    Input rev_fastq     --> 'R2'
To use this script with more pipleines, update the pipelines.json file.
"""

import argparse
import json
import os
import re
import sys

import pandas as pd
import s3fs

def parse_user_input():
    usage = "python create_submission_commands.py   --seed seedfile.txt \
                                                    --s3_output s3://czbiohub-microbiome/My_Lab/My_Name/My_Project/ \
                                                    --commands aegea_launch_commands.sh \
                                                    --pipeline midas \
                                                    --config pipeline.json"

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    # Required
    p.add_argument(
        "-o",
        "--s3_output",
        dest="out_s3path",
        action="store",
        type=str,
        required=True,
        help="S3 Path to the output directory to store the results",
    )
    p.add_argument(
        "-a",
        "--commands",
        dest="aegea_cmd",
        action="store",
        type=str,
        required=True,
        help="name of the file where the aegea commands for each input will be stored",
    )
    p.add_argument(
        "-p",
        "--pipeline",
        dest="pipeline",
        action="store",
        type=str,
        required=True,
        help="name of the pipeline to execute. Should be a valid pipeline, present in the pipelines.json file.",
    )
    p.add_argument(
        "-c",
        "--config",
        dest="config_file",
        action="store",
        type=str,
        default="pipelines.json",
        required=True,
        help="path to the pipelines.json file",
    )

    # Pick one
    p.add_argument(
        "-s",
        "--seed",
        dest="seedfile",
        action="store",
        type=str,
        help='(preferred) comma-delimited file containing 3 columns, named: "sampleName","R1" and "R2"',
    )
    p.add_argument(
        "-i",
        "--s3_input",
        dest="in_s3path",
        action="store",
        type=str,
        help="S3 Path to paired fastq files to be used as inputs in the pipeline",
    )

    # Optional
    ## AWS Logistics
    p.add_argument("--to_aws", dest="xfer", action="store_true")
    p.add_argument("--yes", dest="no_confirmation", action="store_true")
    # p.add_argument("--se", dest="single_ended", action="store_true")
    p.add_argument("--image", dest="image_name", action="store", type=str)
    p.add_argument("--image_version", dest="image_version", action="store", type=str)
    p.add_argument("--queue", dest="queue", action="store", type=str)
    p.add_argument("--script", dest="script", action="store", type=str)
    # p.add_argument(
    #     "--extra",
    #     dest="extraList",
    #     action="append",
    #     type=str,
    #     help="Use this flag to add custom export variables. Can be used multiple times in a single command to add multiple variables. You will be asked for a confirmation (unless --yes specified).",
    # )
    p.add_argument(
        "--force",
        dest="overwrite",
        action="store_true",
        default=False,
        help="Overwrite existing results",
    )

    ## Hardware
    p.add_argument("--memory", dest="memory", action="store", type=int)
    p.add_argument("--core", dest="cpu", action="store", type=int)
    p.add_argument(
        "--storage", dest="storage", action="store", type=int
    )  # the minimum for AWS is 500
    p.add_argument("--retry", dest="max_retries", action="store", type=int)

    ## Tagging
    p.add_argument(
        "-n",
        "--project",
        dest="project",
        action="store",
        type=str,
        help="name of the project that this set of jobs will belong to",
    )
    p.add_argument(
        "-l",
        "--lead",
        dest="team_lead",
        action="store",
        type=str,
        help="name of the team lead / PI responsible for this set of jobs",
    )

    arguments = p.parse_args()
    pipeline = arguments.pipeline.lower()

    with open(arguments.config_file) as pipeline_defaults:
        # default = json.load(pipeline_defaults)
        try:
            params = json.load(pipeline_defaults)[pipeline]
        except KeyError as e:
            print(
            f"[FATAL] This script has not been configured to be used for {pipeline}.\nPlease update {arguments.config_file} and try again."
            )
            sys.exit(1)

    # Ask for required variables
    set_of_params=set()
    if "extras" in params:
        print(f"The {pipeline} pipeline can make use of additional variables, please provide their values below.")
        print(f"NOTE: Not all variables may be required, hit 'return', for no value")
        for variable, description in params["extras"].items():
            user_value = input(f"\t'{variable}' ({description}):\n\t")
            if user_value == "":
                continue
            params.update({variable:user_value})
            set_of_params.add(variable)

    # Set Defaults
    # if default[pipeline]:
    #     data_mount = default[pipeline]["data_mount"]
    #     script_loc = default[pipeline]["script"]
    #     # current setup requires at least 500GB
    #     storage = default[pipeline]["storage"]
    #     memory = default[pipeline]["memory"]
    #     cpu = default[pipeline]["cpu"]
    #     queue = default[pipeline]["queue"]
    #     image = default[pipeline]["image"]
    #     version = default[pipeline]["version"]

    # Default Hard Limits
    min_storage = 500
    max_memory = 488000
    max_cpus = 64

    args = vars(arguments)
    for key, value in params.items():
        if key == "extras":
            continue
        
        set_of_params.add(key)
        if (key in args) and (args[key] is not None):
            if key == "storage":
                params["storage"] = max(min_storage, args["storage"])        
            elif key == "memory":
                params["memory"] = min(max_memory, args["memory"])
            elif key == "cpus":
                params["cpus"] = min(max_cpus, args["cpu"])
            else:
                params.update({key:value})

    if not arguments.no_confirmation:
        while True:
            print("Please confirm that the following variables are set to your satisfaction:")
            print("NOTE: Some values may have been adjusted to satisfy minimum requirements")
            for param in set_of_params:
                print(f"\t{param} = {params[param]}")
            confirmed = input("Satisfied? [Y/n]: ")
            if confirmed == "n":
                print("My apologies, please restart the script and try again... Exiting...")
                sys.exit(0)
            elif confirmed == "Y":
                print("Confirmed. Processing inputs ...")
                break
            else:
                print(
                    f'\n[Invalid Response] "{confirmed}" is not a valid response. Please select either "Y" or "n" (case-sensitive)\nTry again:'
                )

    # remove keys that are in params from args
    for key in set(args) & set(params):
        del args[key]

    args["out_s3path"] = arguments.out_s3path.rstrip("/")

    return args, params

# extra = None
# if len(arguments.extraList) > 0:
#     extraList = map(lambda x: str(x).rstrip(r";$"), arguments.extraList)
#     extra = "".join([f"export {e};" for e in extraList])
#     if not arguments.no_confirmation:
#         while True:
#             confirm = input(
#                 f"""
#         The following statement(s) will be added, as is to the export variables:

#         {extra}

#         Are you sure you wish to continue? [Y/n]
#         """
#             )
#             if confirm == "n":
#                 sys.exit(0)
#             elif confirm == "Y":
#                 print("Confirmed. Processing inputs ...")
#                 break
#             else:
#                 print(
#                     f'\n[Invalid Response] "{confirm}" is not a valid response. Please select either "Y" or "n" (case-sensitive)\nTry again:'
#                 )
# else:
#     extra = ""

def submit_job(
    name,
    fwd,
    submission,
    job,
    rev=None,
):
    job_queue=job["queue"]
    img_version=job["version"]
    docker_image=job["image"]
    job_storage=int(job["storage"])
    job_cpu=int(job["cpu"])
    job_memory=int(job["memory"])
    job_attempts=int(job["retries"])
    job_data_mount=job["data_mount"]
    job_script=job["script"]
    
    s3output=submission["out_s3path"]
    pipeline=submission["pipeline"]
    overwrite=submission["overwrite"]

    aegea_cmd = ""
    memoryPerCore = f"{int(job_memory/(job_cpu * 1000))}G"
    s3output_path = f"{s3output}/{name}"
    if overwrite:
        complete = False
    else:
        complete = fs.exists(f"{s3output_path}/job.complete")

    execute_cmd=""
    if "extras" in job: # Do I need to look for other variables?
        for extra in job["extras"]: # Which variables?
            if extra in job: # Did user provide a value for this variable?
                execute_cmd=f"export {extra}={job[extra]};" # Add the variable to the command

    if not complete:
        execute_cmd = f"{execute_cmd}export coreNum={job_cpu};export memPerCore={memoryPerCore};export fastq1={fwd};export S3OUTPUTPATH={s3output_path}"
        if rev is not None:
            execute_cmd = f"{execute_cmd}; export fastq2={rev}"
        execute_cmd = f"{execute_cmd}; {job_script};"
        aegea_cmd = f"aegea batch submit --retry-attempts {job_attempts} --name {pipeline}__{name} --queue {job_queue} --image {docker_image}:{img_version} --storage {job_data_mount}={job_storage} --memory {job_memory} --vcpus {job_cpu} --command='{execute_cmd}'"

    return aegea_cmd


if __name__ == "__main__":
    fs = s3fs.S3FileSystem(anon=False)
    submission_args, job_args = parse_user_input()

    seed_df = pd.DataFrame()
    if submission_args["seedfile"]:
        try:
            # Comma sep
            seed_df = pd.read_csv(submission_args["seedfile"])
            seed_df["sampleName"]
        except KeyError as e:
            # Tab sep
            seed_df = pd.read_table(submission_args["seedfile"])
            seed_df["sampleName"]
        except KeyError as e:
            print(f"[FATAL] {e}")
            print(seed_df.sample(10))
            sys.exit(1)
        finally:
            print(f"Sample:\n{seed_df.head(5)}")
    elif submission_args["in_s3path"]:
        in_s3path = submission_args["in_s3path"].rstrip("/")
        df = pd.DataFrame()
        df["FilePaths"] = fs.glob(in_s3path + "/*fastq.gz")

        try:
            df["sampleName"] = df["FilePaths"].str.extract(
                r"([a-zA-Z0-9_\-\.]+)_S\d+_L\d+_R[12]_\d+\.fastq\.gz"
            )
            df["Orientation"] = df["FilePaths"].str.extract(
                r"[a-zA-Z0-9_\-\.]+_S\d+_L\d+_(R[12])_\d+\.fastq\.gz"
            )
            df["FilePaths"] = df["FilePaths"].apply(lambda x: "s3://" + x)
            seed_df = df.pivot(
                index="sampleName", columns="Orientation", values="FilePaths"
            )
        except ValueError as e:
            df["sampleName"] = df["FilePaths"].str.extract(
                r"([a-zA-Z0-9_\-\.]+)_S\d+_R[12]_\d+\.fastq\.gz"
            )
            df["Orientation"] = df["FilePaths"].str.extract(
                r"[a-zA-Z0-9_\-\.]+_S\d+_(R[12])_\d+\.fastq\.gz"
            )
            # df["FilePaths"] = df["FilePaths"].apply(lambda x: 's3://' + x)
            seed_df = df.pivot(
                index="sampleName", columns="Orientation", values="FilePaths"
            )
        # except _ as e:
        #     print(f"[FATAL] {e}")
        #     print(df.sample(10))
        #     sys.exit(1)
        finally:
            print(f"Sample:\n{seed_df.head(5)}")

        seed_df.columns.name = None
        seed_df = seed_df.reset_index()
        seedfile = f"{os.path.splitext(submission_args['aegea_cmd'])[0]}.seedfile.csv"
        seed_df.to_csv(seedfile, index=False)
        # seed_df ["commands"] = seed_df.apply(lambda row: submit_job(name = row['sampleName'], fwd = row['R1'], rev = row['R2']), axis=1)
    else:
        print(
            "[FATAL] At least one form of input required. Please either submit a seedfile or a s3 folder with paired fastq.gz files"
        )
        sys.exit(1)

    # Handle special characters in Sample Names
    seed_df["sampleName"] = seed_df["sampleName"].map(lambda x: re.sub(r"\W+", "-", x))
    if "R2" in seed_df.columns:
        commands = seed_df.apply(
            lambda row: submit_job(name=row["sampleName"], fwd=row["R1"], rev=row["R2"], submission=submission_args, job=job_args), axis=1
        ).dropna(axis="index")
    else:
        commands = seed_df.apply(
            lambda row: submit_job(name=row["sampleName"], fwd=row["R1"], submission=submission_args, job=job_args), axis=1
        ).dropna(axis="index")

    commands.to_csv(submission_args['aegea_cmd'], header=False, index=False)
    nrows = commands.shape[0]
    print(f"{nrows} commands written to {submission_args['aegea_cmd']}")

    if submission_args['xfer'] and (nrows > 0):
        aegea_cmd_filename = os.path.basename(submission_args['aegea_cmd'])
        fs.put(submission_args['aegea_cmd'], f"{submission_args['out_s3path']}/{aegea_cmd_filename}")

