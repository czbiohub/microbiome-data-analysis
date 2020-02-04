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
# The $script_name should be aws_ec2_submit_docker_basecall.txt
#
# bash aws_ec2_launch_gpu_instance $iam_role Path_to_$script_name
#
# This script outputs a json file in the current directory.
#
# Revision History
# 2019.09.08 Created by Brian Yu
#####################################################################

# Take iam_role as the first argument into the script
iam_role=$1
# script file name is the second argument. Script name must end in .txt
script_name=$2

# Define parameters used for starting the EC2 instance
pem_name=czbiohub_microbiome # This is what you use only if you want to ssh into the instance
nanopore_ami=ami-0ddba16a97b1dcda5 # This is a public deep learning ami
new_block_device=DeviceName=/dev/sda1,Ebs={VolumeSize=1000} # Adds storage at root so that the docker image can be installed

# Information for tracking purposes.
# To ssh into the instance manually - aegea ssh ubuntu@$instance_name -i $pem_name.pem
# Although you must have czbiohub_microbiome.pem on your path and do - chmod 400 czbiohub_microbiome.pem
instance_name=Nanopore_Basecall_$( date +"%Y%m%d_%H%M%S" )
new_tags=ResourceType=instance,Tags=[{Key=Name,Value=$instance_name},{Key=Project,Value=Microbiome_Nanopore}]

# Start Instances. --user-data must end in .txt and contain #!/bin/bash in first row
# shutdown behavior is stop or terminate. In this case use terminate.
# In the script use "shutdown -h now", which will terminate instance p3.2xlarge

echo "COMMAND TO RUN ::: aws ec2 run-instances --image-id $nanopore_ami --instance-type p3.2xlarge --count 1 \
--instance-initiated-shutdown-behavior "terminate" --block-device-mapping $new_block_device \
--key-name $pem_name --iam-instance-profile Name=$iam_role \
--user-data file://$script_name --tag-specifications $new_tags"

aws ec2 run-instances --image-id $nanopore_ami --instance-type p3.2xlarge --count 1 \
--instance-initiated-shutdown-behavior "terminate" --block-device-mapping $new_block_device \
--key-name $pem_name --iam-instance-profile Name=$iam_role \
--user-data file://$script_name --tag-specifications $new_tags > $instance_name.json

# Update instance information with describe-Instances
sleep 60
instance_id=$( grep InstanceId $instance_name.json | cut -d\" -f4 )
aws ec2 describe-instances --instance-id $instance_id > $instance_name.json

# The problem with aegea launch is it does not allow user defined device name to be /dev/sda1
# The function that handles the --storage option is the get_bdm function in aegea/aegea/util/aws/__init__.py
