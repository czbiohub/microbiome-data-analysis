#!/usr/bin/env python3

import boto3
import pandas as pd
import time
import sys
import logging

def get_jobs_list(batch,queue_name):
    lodf=list()
    response = batch.list_jobs(jobQueue=queue_name)
    lodf.append(response_to_df(response))
    while "nextToken" in response:
        response = batch.list_jobs(jobQueue=queue_name, nextToken=response["nextToken"])
        lodf.append(response_to_df(response))

    return pd.concat(lodf)

def response_to_df(response):
    return pd.DataFrame(response['jobSummaryList'])

def runtime_in_hours(millis, now):
    millis = int(now) - int(millis)
    # seconds=(millis/1000)%60
    # seconds = int(seconds)
    # minutes=(millis/(1000*60))%60
    # minutes = int(minutes)
    hours=(millis/(1000*60*60))%24
    return hours

def terminate_selected_jobs(batch, df):
    return df["jobId"].apply(lambda x: batch.terminate_job(jobId=x,reason="User Initiatied Timeout Termination"))

def main(queue,max_runtime_hrs, output):
    """[summary]
    main("sonnenburg__spot100",5, "TEST.terminate_jobs.csv")

    Args:
        queue ([type]): [description]
        max_runtime_hrs ([type]): [description]
        output ([type]): [description]
    """
    batch = boto3.client('batch')
    current_time = time.time_ns() // 1000000 

    df = get_jobs_list(batch, queue)
    njobs, _ = df.shape
    if njobs == 0:
        logging.info(f"No currently running jobs found in queue ({queue}). Exiting.")
        return

    df["runtime_hrs"] = df["startedAt"].apply(lambda x: runtime_in_hours(x, current_time))
    selected_df = df.query(f"runtime_hrs > {max_runtime_hrs}")[["jobId","jobName","runtime_hrs"]]
    term_njobs, _ = selected_df.shape
    
    if term_njobs == 0:
        logging.info(f"No jobs running beyond the maximum runtime hours ({max_runtime_hrs}) found. Exiting.")
        return
    
    while True:
        logging.info(f"Found {njobs} total jobs running, {term_njobs} jobs can be terminated beacuse they exceed your time out criteria of {max_runtime_hrs} hrs.")
        logging.warning("*WARNING*: This action CAN NOT be undone. You have been warned!")
        terminate_jobs = input("Would you like to terminate these jobs?[Y/n]: ")
        if terminate_jobs == "Y":
            terminate_selected_jobs(batch, selected_df)
            break
        elif terminate_jobs == "n":
            logging.info(f"You selected [{terminate_jobs}]")
            logging.info(f"No actions will be performed.")
            break
        else:
            logging.info(f"You selected [{terminate_jobs}]. This is not a valid response. Please choose 'Y' to terminate jobs, or 'n' to just save the list of jobs. (The response validation is case-sensitive)")

    selected_df.to_csv(output, index=False)
    return 

if __name__=="__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )
    queue=sys.argv[1]
    max_runtime_hrs=sys.argv[2]
    output=sys.argv[3]
    main(queue, max_runtime_hrs, output)