# Nanopore GPU Basecall, Demux and Trim

The `nanopore_basecall_demux_gpu` docker container can be used to basecall, demultiplex, and trim adapters for Nanopore MinION runs. It outputs one `.fastq.gz` file per library if a run is barcoded. Currently, the script has only been tested with `MIN106` flow cells and Native Ligation Barcode expansion kits `EXP-NBD104` and/or `EXP-NBD114`. The docker image should be run on AWS EC2 and the output files will be saved to AWS S3. `Guppy v4.0` is used to perform basecall and demux. `Porechop` is used to perform adapter trimming.

For questions, please contact [Brian Yu](brian.yu@czbiohub.org)

## Using `nanopore_basecall_gpu` docker image

Upload the Nanopore sample sheet to the appropriate AWS S3 location
```
aws s3 cp ${path_to_local_sample_sheet.csv} ${s3_sample_sheet_path}
```
Make a copy of the `aws_ec2_launch_gpu_instance_template.sh` file to your local directory, edit `line 25` team leader value `Value=FirstName_LastName` to reflect your team leader.

Run the docker image
```
bash aws_ec2_launch_gpu_instance_template.sh $iam_role $nanopore_run_name $s3_path_to_run_folder
```

The script should spend ~60 seconds and produce an output file containing information about the instance in `json` format.

The docker image must use an AWS EC2 instance of the `p3.2xlarge` type, with `--block-device-mapping=DeviceName=/dev/sda1,Ebs={VolumeSize=1000}` and this public deep learning AMI `ami-0ddba16a97b1dcda5`.


## Building `nanopore_basecall_gpu` from command line

The current docker image is built to use a specific AWS S3 location for sample sheet storage. The sample sheet template is included `samplesheet_template.csv`. If you do not have the proper credentials and/or want to use this docker image with another S3 location, you will need rebuild the docker image. To do so, edit `nanopore_basecall_demux.sh` `line 78` to reflect the new sample sheet location.
```
python sample_sheet_process.py -t ${coreNum} -i ${RUN_NAME}.csv -p ${s3_sample_sheet_path}
```

You will also need to edit `Makefile` `line 3` to correspond to your docker account.
```
REGISTRY=$YOUR_DOCKER_ACCOUNT_NAME/$(NAME)
```
Make sure you are connected to your docker account, then run the `Makefile`
```
make all
```

## Files used by `nanopore_basecall_gpu` docker image

`Dockerfile` describes what is installed in the docker image.

`Makefile` is used to build and update the docker image.

`install_packages.R` includes R packages installed in the docker image. These R packages are used by the MinIONQC.R script from [here](https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R).

`LATEST` contains the name of the latest version of the `nanopore_basecall_gpu` docker image.

`nanopore_basecall_demux.sh` is the script that is used to perform the actual basecalling and demultiplexing functionalities. It calls `sample_sheet_process.py` to perform adapter trimming using [porechop](https://github.com/rrwick/Porechop).

`aws_ec2_launch_gpu_instance_template.sh` is a bash script that can be used to programmatically launch the proper EC2 instance and run the docker images for a particular MinION Nanopore run.

`samplesheet_template.csv` is an example sample sheet that is compatible with the docker image.
