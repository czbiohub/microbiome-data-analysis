#!/bin/bash

####################################################################
# aws_ec2_launch_gpu_instance.sh
#
# Run this script as the entry point to nanopore demux using Bryan Merrill's
# docker image bmerrill9/nanopore_basecall_demux_gpu:latest
#
# It uses the pem file czbiohub_microbiome.pem
# A public deep learning ami ami-0ddba16a97b1dcda5
# It can be invoked as below.
# You need to know an IAM role to use
#####################################################################

iam_role=$1
RUN_FOLDER_NAME=$2
s3_PATH=$3

# Define parameters used for starting the EC2 instance
nanopore_ami=ami-0ddba16a97b1dcda5 # This is a public deep learning ami
new_block_device=DeviceName=/dev/sda1,Ebs={VolumeSize=1000} # Adds storage at root so that the docker image can be installed

# Information for tracking purposes.
instance_name=Nanopore_Basecall_$( date +"%Y%m%d_%H%M%S" )
new_tags=ResourceType=instance,Tags=[{Key=Name,Value=$instance_name},{Key=Project,Value=Nanopore_Basecall_Demux},{Key=Team_Leader,Value=FirstName_LastName}]

# Check to see if the sample sheet exists before launching the instance
exist_flag=$( aws s3 ls s3://czb-seqbot/nanopore/nanopore-samplesheets/${RUN_FOLDER_NAME}.csv )
if [ -z "${exist_flag}" ]
then
  echo "Sample sheet ${RUN_FOLDER_NAME}.csv is not found in s3://czb-seqbot/nanopore/nanopore-samplesheets/"
  exit
else

  # Create a file to use as --user-data
  echo \#\!/bin/bash > user_data_file.txt
  echo -e "docker container run --runtime=nvidia -e NANOPORE_RUNPATH=${s3_PATH} -e RUN_NAME=${RUN_FOLDER_NAME} --name=nanopore_basecall_demux brianyu2010/nanopore_basecall_demux_gpu:latest \"./nanopore_basecall_demux.sh\"" >> user_data_file.txt
  echo -e "docker wait nanopore_basecall_demux" >> user_data_file.txt
  echo -e "shutdown -h now" >> user_data_file.txt

  echo "COMMAND TO RUN ::: aws ec2 run-instances --image-id $nanopore_ami --instance-type p3.2xlarge --count 1 \
  --instance-initiated-shutdown-behavior ""terminate"" --block-device-mapping $new_block_device \
  --iam-instance-profile Name=$iam_role \
  --user-data file://user_data_file.txt --tag-specifications $new_tags"

  aws ec2 run-instances --image-id $nanopore_ami --instance-type p3.2xlarge --count 1 \
  --instance-initiated-shutdown-behavior terminate --block-device-mapping $new_block_device \
  --iam-instance-profile Name=$iam_role \
  --user-data file://user_data_file.txt --tag-specifications $new_tags > $instance_name.json

  # file://$script_name
  # Update instance information with describe-Instances
  sleep 60
  instance_id=$( grep InstanceId $instance_name.json | cut -d\" -f4 )
  aws ec2 describe-instances --instance-id $instance_id > $instance_name.json

fi
