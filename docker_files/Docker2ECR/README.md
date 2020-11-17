# Docker --> ECR

## Why?

Needed to move my docker images to ECR, so I could continue using my current AWS Batch based infrastructure with minimal downtime.

## How does it work?

The script uses `boto3` and `docker` python packages. The script first searches your local enviroment for the image that you wish to move, if it can't find it, it'll search public repos (based on your input). Once an image is found, it'll use your AWS credentials to upload a copy to your ECR.

Special care has been taken to add images to specific namespaces. In it's current form, the script supports two Stanford labs and requires that you be a part of one of these to use it (although, it is easy enough to bypass/update this behaviour).

## Usage

python docker2ecr.py local_image_name ecr_image_name image_tag lab_of

### Example 1

When you don't have an image on your system. Script downloads the docker `image:tag` based on input

```{bash}
   python docker2ecr.py sunitjain/midas midas 20191119142538 fischbach
```

### Example 2

When you have an image on your system. Script uses local docker `image:tag` based on input

```{bash}
   python docker2ecr.py midas midas 20191119142538 fischbach
```

## Requirements

- Make sure you have configured access to AWS. The script uses your default profile and `us-west-2` region.
- You will also need to log in to docker via your terminal. You may download and install docker desktop and type `docker login` in your terminal and follow the prompts.
- Python packages:
  - boto3
  - docker

## Credits

Borrowed from : https://raw.githubusercontent.com/AlexIoannides/py-docker-aws-example-project/master/deploy_to_aws.py