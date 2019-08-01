# FiBo

This repo is a frozen dockerized fork of the original repo:

- [FitnessBrowser](https://github.com/pflu-sbw25/FitnessBrowser).
- Branch: `master`
- commit: `9c85d32a4f0b9baca7e27604d20cbbc0a45bce33`
- dated: `March 21, 2019`

Currently, it supports functionality that is rather limited to my projects, but, additional support may be added at a later date depending on the interest.

## MultiCodes.pl

This script has been wrapped into `run_multiCodes.sh` and expects at least an S3 path to a compressed fastq (with the extension `.fastq.gz`) and an S3 output location. You may also provide additional values such as the `MIN_QUALITY`, `INDEX_DIR` and `INDEX_FILE`.

Submission using `aegea`:

```bash
aegea batch submit \
    --watch --retry-attempts 1 \
    --queue your-favorite-queue \
    --image sunitjain/fibo:latest \
    --storage /mnt=500 --memory 4000 --vcpus 2 \
    --command="export MIN_QUALITY=20;export S3OUTPUTPATH=s3://some_bucket/FiBo/20190801_BarSeq/Sample_Prefix; export S3INPUTPATH=s3://some_bucket/FiBo/20190801_BarSeq/FastQ/Sample_Prefix.fastq.gz;./run_multiCodes.sh"
```

## Questions/Beers

[Sunit Jain](microbiome.ninja) (dev)