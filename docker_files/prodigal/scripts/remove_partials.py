import sys
import argparse
from Bio import SeqIO

infile = sys.argv[1]
prefix = sys.argv[2]
outdir = sys.argv[3].rstrip("/")
suffix = sys.argv[4]

outfile = f"{outdir}/{prefix}.{suffix}"

complete = 0
incomplete = 0
with open(outfile, "w") as output_handle:
    for record in SeqIO.parse(infile, "fasta"):
        if "partial=00" in record.description:
            complete += 1
            record.id = f"{prefix}__{record.id}"
            SeqIO.write(record, output_handle, "fasta")
        else:
            incomplete += 1
            # print(f'{record.id} {record.description}')

print(f"{prefix}\t{complete}\t{incomplete}")
