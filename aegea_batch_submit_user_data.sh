#!/bin/bash -c 

for i in "$@"; do 
    eval "$i"
done 

aegea.batch 

set -a 

if [ -f /etc/environment ]; then 
    source /etc/environment
fi 
 
if [ -f /etc/default/locale ];then 
    source /etc/default/locale
fi 

set +a 

if [ -f /etc/profile ];then 
    source /etc/profile
fi 

set -euo pipefail 

apt-get update -qq apt-get install -qqy --no-install-suggests --no-install-recommends httpie awscli jq psmisc lsof 

iid=$(http http://169.254.169.254/latest/dynamic/instance-identity/document) 
aws configure set default.region $(echo "$iid" | jq -r .region) 

az=$(echo "$iid" | jq -r .availabilityZone) 

echo "Creating volume" >&2 
vid=$(aws ec2 create-volume --availability-zone $az --size 500 --volume-type st1 | jq -r .VolumeId) 

aws ec2 create-tags --resource $vid --tags Key=aegea_batch_job,Value=$AWS_BATCH_JOB_ID 

echo "Setting up SIGEXIT handler" >&2 
trap "cd / 
 fuser /mnt >&2 || echo \"Fuser exit code \$?\" >&2
 lsof /mnt | grep -iv lsof | awk '{print \$2}' | grep -v PID | xargs kill -9 || echo \"LSOF exit code \$?\" >&2
 sleep 3

 umount /mnt || umount -l /mnt
 aws ec2 detach-volume --volume-id $vid

 let try=1
 sleep 10

 while ! aws ec2 describe-volumes --volume-ids $vid | jq -re .Volumes[0].Attachments==[]; do 
    if [[ \$try -gt 2 ]]; then 
        echo \"Forcefully detaching volume $vid\" >&2
        aws ec2 detach-volume --force --volume-id $vid
        sleep 10
        echo Deleting volume $vid >&2
        aws ec2 delete-volume --volume-id $vid
        exit
    fi
    sleep 10
    let try=\$try+1
 done
 
 echo \"Deleting volume $vid\" >&2
 
 aws ec2 delete-volume --volume-id $vid" EXIT 
 
echo "Waiting for volume $vid to be created" >&2 
 
while [[ $(aws ec2 describe-volumes --volume-ids $vid | jq -r .Volumes[0].State) != available ]]; do 
   sleep 1
done 
 
let pid=$$ 
 
echo "Finding unused devnode for volume $vid" >&2 

let delay=2+$pid%5 

for try in {1..100}; do 
    if [[ $try == 100 ]]; then 
        echo "Unable to mount $vid on $devnode"
        exit 1
    fi

    if [[ $try -gt 1 ]]; then 
        sleep $delay
    fi
 
    devices=$(echo /dev/xvd* /dev/xvd{a..z} /dev/xvd{a..z} | sed 's= =\n=g' | sort | uniq -c | sort -n | grep ' 2 ' | awk '{print $2}')
    let devcnt=${#devices}/10+1
    let ind=$pid%devcnt
    devnode=${devices:10*$ind:9}
    aws ec2 attach-volume --instance-id $(echo "$iid" | jq -r .instanceId) --volume-id $vid --device $devnode || continue
    break
done 
 
echo "Waiting for volume $vid to attach on $devnode" >&2 
let delay=5+$pid%11 
sleep $delay 

let try=1 
let max_tries=32 

while [[ $(aws ec2 describe-volumes --volume-ids $vid | jq -r .Volumes[0].State) != in-use ]]; do 
    if [[ $try == $max_tries ]]; then 
        break
    fi
    
    let foo=1+$try%5
    let delay=2**$foo+$pid%11
    sleep $delay
done 
 
while [[ ! -e $devnode ]]; do 
    sleep 1
done 
 
echo "Making filesystem on $devnode" >&2 
mkfs.ext4 $devnode 

echo "Mounting $devnode" >& 2 

mount $devnode /mnt 

echo "Devnode $devnode mounted" >& 2 

./run_sourmash.sh