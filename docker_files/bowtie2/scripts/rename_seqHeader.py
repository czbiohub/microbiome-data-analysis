import sys
import argparse
from Bio import SeqIO

infile = sys.argv[1]
prefix = sys.argv[2]
outdir = sys.argv[3].rstrip("/")
suffix = sys.argv[4]

outfile = f"{outdir}/{prefix}.{suffix}"

complete = incomplete = 0
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

# python3 scripts/01_filter_alignments.py \
#     -bam /mnt/data/PacBio/minimap/Mixture/alignment/bam/m64069_200107_192940.subreads.bam.ccs_vs_SCv2_index.sortedByCoord.bam \
#     -fasta /mnt/data/PacBio/minimap/Mixture/db/20190911_scv2.fna \
#     -outdir /mnt/data/PacBio/minimap/Mixture/alignment/NinjaMapPB_0_0/ \
#     -bin /mnt/data/PacBio/minimap/Mixture/db/20190911_scv2.bin.tsv \
#     -se \
#     -min_pid 0 \
#     -min_paln 0 &> /mnt/data/PacBio/minimap/Mixture/alignment/NinjaMapPB_0_0/01_filter_alignments.log
