#!/usr/bin/env python3
"""
docker2ecr.py
~~~~~~~~~~~~~~~~
Borrowed from : https://raw.githubusercontent.com/AlexIoannides/py-docker-aws-example-project/master/deploy_to_aws.py

A simple script that deploys an existing local or remote docker image
repository to your ECR account.

USAGE:
python docker2ecr.py local_image_name ecr_image_name image_tag fischbach lab_of

Example 1: When you don't have an image on your system. Script downloads
the docker image:tag based on input

   $ python docker2ecr.py sunitjain/midas midas 20191119142538 fischbach

Example 2: When you have an image on your system. Script uses local docker image:tag 
based on input

   $ python docker2ecr.py midas midas 20191119142538 fischbach
"""

import base64
import logging
import os
import sys

import boto3
import docker


def main(local_image_name, ecr_image_name, image_tag, lab_of):
    """Get prebuilt Docker image, push to AWS and update ECS service.

    Args:
        local_image_name (string): image name (local or on dockerhub or elsewhere; no tags)
        ecr_image_name (string): name of image on ECR (without the account prefixes; see examples)
        image_tag (string): image tag to copy into ECR
        lab_of (string): Which lab does this image belong to? Only "fischbach" and "sonnenburg" are approved.

    Returns:
        None
    """
    logging.info("Making sure all tags and permissions are in order ...")
    ecr_repo_name, repo_tags = namespace_check(lab_of, ecr_image_name)

    docker_client = docker.from_env()
    image, cleanup_required = get_image(docker_client, local_image_name, image_tag)

    # get AWS ECR login token
    logging.info("Getting ECR Authorization Token ...")
    ecr_client = boto3.client("ecr", region_name="us-west-2")

    # Create or Get ECR repo details
    try:
        ecr_repo = ecr_client.create_repository(
            repositoryName=ecr_repo_name, tags=repo_tags
        )["repository"]
    except ecr_client.exceptions.RepositoryAlreadyExistsException as e:
        ecr_repo = ecr_client.describe_repositories(repositoryNames=[ecr_repo_name])[
            "repositories"
        ][0]

    ecr_credentials = ecr_client.get_authorization_token()["authorizationData"][0]
    ecr_username = "AWS"
    ecr_password = (
        base64.b64decode(ecr_credentials["authorizationToken"])
        .replace(b"AWS:", b"")
        .decode("utf-8")
    )
    ecr_url = ecr_credentials["proxyEndpoint"]

    # tag image for AWS ECR
    ecr_uri = ecr_repo["repositoryUri"]
    image.tag(ecr_uri, tag=image_tag)

    logging.info("Pushing image to ECR ...")
    # get Docker to login/authenticate with ECR
    docker_client.login(username=ecr_username, password=ecr_password, registry=ecr_url)
    # push image to AWS ECR
    for line in docker_client.images.push(
        repository=ecr_uri,
        tag=image_tag,
        auth_config={"username": ecr_username, "password": ecr_password},
        stream=True,
        decode=True,
    ):
        logging.info(line)

    if cleanup_required:
        docker_client.images.remove(image=f"{ecr_uri}:{image_tag}")
        docker_client.images.remove(image=f"{local_image_name}:{image_tag}")

    logging.info(f"You may try: '{ecr_uri}:{image_tag}'")

    return None


def namespace_check(lab_of, ecr_image_name):
    approved_team_leaders = {
        "fischbach": "Michael Fischbach",
        "sonnenburg": "Justin Sonnenburg",
    }

    ## Team Leader spell check
    assert (
        lab_of.lower() in approved_team_leaders
    ), f"Your lab '{lab_of}' is not approved for an ECR repo. You MUST belong to one of the approved labs. Please contact sunit at czbiohub."

    team_leader = approved_team_leaders[lab_of.lower()]

    approved_namespaces = {
        "Michael Fischbach": "fischbach_lab",
        "Justin Sonnenburg": "sonnenburg_lab",
    }

    repo_namespace = f"{approved_namespaces[team_leader]}/{ecr_image_name}"
    repo_tags = generate_repo_tags(team_leader)

    return repo_namespace, repo_tags


def generate_repo_tags(team_leader):
    ownership = {
        "Project": "Various - Microbiome Initiative",
        "Team Leader": f"{team_leader}",
    }
    repo_tags = list()
    for key, value in ownership.items():
        repo_tags.append({"Key": f"{key}", "Value": f"{value}"})

    return repo_tags


def get_image(client, repo_name, tag="latest"):
    """Acquire the specified image object, locally or via the internets

    Args:
        client (obj): docker client object
        repo_name (string): name of the image
        tag (optional, string): specific tag (default: "latest")

    Returns:
        [type]: [description]
    """
    image = f"{repo_name}:{tag}"
    cleanup = False
    try:
        # Image is already built and exists locally
        logging.info(f"Searching for image '{image}' locally ...")
        image = client.images.get(image)
    except docker.errors.ImageNotFound as notfound:
        # Image is already built and exists in a repo
        logging.info(
            f"Could not find image '{image}' locally. Will attempt to download from the internets ..."
        )
        image = client.images.pull(repo_name, tag=tag)
        cleanup = True
    except docker.errors.APIError as apierror:
        logging.error(
            f"Could not find image '{image}' locally or on the internets. Are you sure it exists?"
        )
        print(apierror)
        sys.exit(1)

    return image, cleanup


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    ## Example of a local docker image
    # local_image_name = "readqc"
    # ecr_image_name = "readqc"
    # image_tag = "20200911001353"
    # lab_of = "fischbach" # or "sonnenburg" for Sonnenburg Lab

    ## Example of a remote docker image
    # local_image_name = "sunitjain/midas"
    # ecr_image_name = "midas"
    # image_tag = "20191119142538"
    # lab_of = "fischbach" # or "sonnenburg" for Sonnenburg Lab

    local_image_name = sys.argv[1]
    ecr_image_name = sys.argv[2]
    image_tag = sys.argv[3]
    # lab_of: must be one of "fischbach" or "sonnenburg"
    lab_of = sys.argv[4]

    main(local_image_name, ecr_image_name, image_tag, lab_of)
