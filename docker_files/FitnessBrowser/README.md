# FiBo

## Contents

This repo is a frozen dockerized fork of the original repo:

- [FitnessBrowser](https://github.com/pflu-sbw25/FitnessBrowser).
- Branch: `master`
- commit: `9c85d32a4f0b9baca7e27604d20cbbc0a45bce33`
- dated: `March 21, 2019`

Currently, it supports functionality that is rather limited to my projects, but, additional support may be added at a later date depending on the interests.

### Changes made to the original repos

`submitter.pl` and `BarSeqTest.pl` scripts were modified by changing their redirect operators from `>& logfile` to a more explicit and compatible `> logfile 2>&1`.

## Before Using Dockerfile

The tool has been installed in `/mnt/FitnessBrowser`, so, in order to call the scripts, you'll need to provide the full path to the scipt, for example: `/mnt/FitnessBrowser/bin/MultiCodes.pl`.

The default workdir for the container is `/data`. This is where you start once the container launches.

## Packaged pipeline

Alternatively, if you use AWS, feel free to use the `bash` script(s) to run specific jobs. If you can't find a wrapper associated with your use case, please feel free to submit an issue or a pull request.

### Use Case: Running a job by sets of samples

Assumption: all the samples were run on the same plate.

This script has been wrapped into `run_feba.sh` and expects at least an S3 path to a compressed fastq (with the extension `.fastq.gz`) and an S3 output location. You may also provide additional values such as the `MIN_QUALITY`, `INDEX_DIR` and `INDEX_FILE`.

Submission using `aegea`:

Use the `create_submission_commands.py` to generate `aegea` submission commands for multiple sets. The script expects a tab-delimited file with the following columns (with the same headers):

- **Set_Name** : An experimental grouping of samples that need to be compared.
- **Sample_Name** : The first `Sample_Name` in each set will be treated as `Time0` or `baseline`.
- **Index_Name** : Two samples in the same set **MAY NOT** have the same `Index_Name`. This version only supports index names in the format IT`XXX` where, `XXX` is a `0` padded number from `001` to `096`.
- **Fastq_S3_Path** : S3 path to demux-ed sequence file. Is expected to be compressed with extension `.fastq.gz`

```bash
aegea batch submit \
    --watch --retry-attempts 1 \
    --queue your-favorite-queue \
    --image sunitjain/fibo:latest \
    --storage /mnt=500 --memory 4000 --vcpus 2 \
    --command="export MIN_QUALITY=20;export S3OUTPUTPATH=s3://some_bucket/FiBo/20190801_BarSeq/Sample_Prefix; export S3INPUTPATH=s3://some_bucket/FiBo/20190801_BarSeq/FastQ/Sample_Prefix.fastq.gz;./run_multiCodes.sh"
```

## Packaged Pipeline

```bash
# Setup Organism Folder
aws s3 sync s3://czbiohub-microbiome/ReferenceDBs/BarSeq/FitnessBrowser/Organisms/bTheta g/bTheta
# Select samples to be processed together; Get this from the user.
ln -s ../Experiments/Plate1_Com_T12.txt .
cat ../SampleSheets/19*index.csv | grep -w -f Plate1_Com_T12.txt | sort > setup.txt
# Download and Organize fastq files
cat setup.txt | parallel --colsep ',' "aws s3 cp {3} fastq_files/{2}_Plate1_Com_T12.fastq.gz" &> Plate1_Com_T12.download.log &
# Create a index line
INDICES=$(cut -d , -f 2 setup.txt | tr "\n" ":" | sed "s/:$//")
# Create a description line (replace T0 sample with the keyword 'Time0')
cut -d , -f 1 setup.txt | tr "\n" ":" | sed "s/:$//"; echo ''
# RUN JOB
perl /mnt/FitnessBrowser/bin/BarSeqTest.pl -bs3 -org bTheta -fastqdir fastq_files -index IT077:IT089:IT083:IT053:IT059:IT004:IT010:IT016:IT022:IT028:IT034:IT040:IT046:IT052:IT058:IT064:IT070:IT076:IT082:IT088:IT094:IT047:IT065:IT071 -desc Time0:Time0:Time0:100_Btmutpool_rep1_T12:100_Btmutpool_rep2_T12:Com1_100_Btmutpool_rep1_T12:Com1_100_Btmutpool_rep2_T12:Com2_100_Btmutpool_rep1_T12:Com2_100_Btmutpool_rep2_T12:Com3_100_Btmutpool_rep1_T12:Com3_100_Btmutpool_rep2_T12:Com4_100_Btmutpool_rep1_T12:Com4_100_Btmutpool_rep2_T12:Com5_100_Btmutpool_rep1_T12:Com5_100_Btmutpool_rep2_T12:Com6_100_Btmutpool_rep1_T12:Com6_100_Btmutpool_rep2_T12:Com7_100_Btmutpool_rep1_T12:Com7_100_Btmutpool_rep2_T12:Com8_100_Btmutpool_rep1_T12:Com8_100_Btmutpool_rep2_T12:Com9__100_Btmutpool_rep2_T12:PBS_rep1_T12:PBS_rep2_T12 &> Plate1_Com_T12.run.log
# Compress the results
tar cvzf Plate1_Com_T12.tar.gz g
# Ship the results to S3
aws s3 cp Plate1_Com_T12.tar.gz s3://czbiohub-microbiome/Fischbach_Lab/Mikhail_Iakiviak/20190824_BarSeq/Results/
```

## Questions/Beers

[Sunit Jain](microbiome.ninja) (dev)

```bash
# Setup Organism Folder
aws s3 sync s3://czbiohub-microbiome/ReferenceDBs/BarSeq/FitnessBrowser/Organisms/bTheta g/bTheta
mkdir -p g/bTheta/results

# Select samples to be processed together; Get this from the user.
ln -s ../Experiments/Plate1_Com_T12.txt .
cat ../SampleSheets/19*index.csv | grep -w -f Plate1_Com_T12.txt | sort > setup.txt
# Download and Organize fastq files
cat setup.txt | parallel --colsep ',' "aws s3 cp {3} fastq_files/{1}__{2}.fastq.gz" &> Plate1_Com_T12.download.log &
# for each fastq:
cat setup.txt | parallel --colsep ',' "zcat fastq_files/{1}__{2}.fastq.gz | /mnt/FitnessBrowser/bin/MultiCodes.pl -bs3 -minQuality 0 -index {2} -out g/bTheta/results/{1}__{2}" &> multicodes.log
# Once all have finished:
# Edit the codes files with the updated the index numbers
CODE_FILES=$(find . -name '*.codes' | tr "\n" " ")
/mnt/FitnessBrowser/bin/combineBarSeq.pl g/bTheta/results/Plate123_BFD_10X g/bTheta/pool.n10 ${CODE_FILES}

# Create the exps file
    # See format below.

/mnt/FitnessBrowser/bin/BarSeqR.pl -org bTheta -exps g/bTheta/barseqtest/exps_table  -pool g/bTheta/pool.n10 -indir g/bTheta/barseqtest -outdir g/bTheta/barseqtest -genes g/bTheta/genes.GC test

```

```bash
# cat ../Plate3_Com_T48/g/bTheta/barseqtest/exps_table
SetName	Index	Description	Date_pool_expt_started	Sample_Name	Location
Plate123_BFD_50X	IT101	Time0	somedate	5X_T0	s3://czb-seqbot/fastqs/190721_NB501938_0142_AHFFCYBGXB/5X_T0_S41_R1_001.fastq.gz
Plate123_BFD_50X	IT100	Time0	somedate	Btmutpool_T0	s3://czb-seqbot/fastqs/190721_NB501938_0142_AHFFCYBGXB/Btmutpool_T0_S44_R1_001.fastq.gz
Plate123_BFD_50X	IT054	50X_repA_T12H	somedate	50X_repA_T12H	s3://czb-seqbot/fastqs/190721_NB501938_0142_AHFFCYBGXB/50X_repA_T12H_S45_R1_001.fastq.gz
Plate123_BFD_50X	IT006	50X_repA_T24H	somedate	50X_repA_T24H	s3://czb-seqbot/fastqs/190721_NB501938_0143_AHMMCLBGXB/50X_repA_T24H_S41_R1_001.fastq.gz
Plate123_BFD_50X	IT060	50X_repA_T48H	somedate	50X_repA_T48H	s3://czb-seqbot/fastqs/190721_NB501938_0143_AHMMCLBGXB/50X_repA_T48H_S93_R1_001.fastq.gz
Plate123_BFD_50X	IT066	50X_repB_T12H	somedate	50X_repB_T12H	s3://czb-seqbot/fastqs/190721_NB501938_0142_AHFFCYBGXB/50X_repB_T12H_S46_R1_001.fastq.gz
Plate123_BFD_50X	IT018	50X_repB_T24H	somedate	50X_repB_T24H	s3://czb-seqbot/fastqs/190721_NB501938_0143_AHMMCLBGXB/50X_repB_T24H_S42_R1_001.fastq.gz
Plate123_BFD_50X	IT072	50X_repB_T48H	somedate	50X_repB_T48H	s3://czb-seqbot/fastqs/190721_NB501938_0143_AHMMCLBGXB/50X_repB_T48H_S94_R1_001.fastq.gz
Plate123_BFD_50X	IT078	50X_repC_T12H	somedate	50X_repC_T12H	s3://czb-seqbot/fastqs/190721_NB501938_0142_AHFFCYBGXB/50X_repC_T12H_S47_R1_001.fastq.gz
Plate123_BFD_50X	IT030	50X_repC_T24H	somedate	50X_repC_T24H	s3://czb-seqbot/fastqs/190721_NB501938_0143_AHMMCLBGXB/50X_repC_T24H_S43_R1_001.fastq.gz
Plate123_BFD_50X	IT084	50X_repC_T48H	somedate	50X_repC_T48H	s3://czb-seqbot/fastqs/190721_NB501938_0143_AHMMCLBGXB/50X_repC_T48H_S95_R1_001.fastq.gz
Plate123_BFD_50X	IT090	50X_repD_T12H	somedate	50X_repD_T12H	s3://czb-seqbot/fastqs/190721_NB501938_0142_AHFFCYBGXB/50X_repD_T12H_S48_R1_001.fastq.gz
Plate123_BFD_50X	IT042	50X_repD_T24H	somedate	50X_repD_T24H	s3://czb-seqbot/fastqs/190721_NB501938_0143_AHMMCLBGXB/50X_repD_T24H_S44_R1_001.fastq.gz
Plate123_BFD_50X	IT096	50X_repD_T48H	somedate	50X_repD_T48H	s3://czb-seqbot/fastqs/190721_NB501938_0143_AHMMCLBGXB/50X_repD_T48H_S96_R1_001.fastq.gz
```
