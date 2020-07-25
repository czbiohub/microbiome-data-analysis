#!/usr/bin/env python3
import sys
import argparse
from Bio import SeqIO

infile = sys.argv[1]
filter_list = sys.argv[2]
outfile = sys.argv[3]

# print(f'#{outfile}')
filter_index = set()
with open(filter_list,"r") as filterFile:
    for line in filterFile:
        filter_index.add(line.strip())
# print(filter_index)

matches = 0
with open(outfile, "w") as output_handle:
    for record in SeqIO.parse(infile, "fasta"):
        if record.id in filter_index:
            matches += 1
            SeqIO.write(record, output_handle, "fasta")

print(f'Found {matches} matches out of {len(filter_index)} items in list')