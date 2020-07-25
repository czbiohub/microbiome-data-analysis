#!/usr/bin/env python3
import sys
from Bio import SeqIO

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


if __name__ == "__main__":
    fastafile = sys.argv[1]
    batch_size = int(sys.argv[2])
    prefix = sys.argv[3]
    print(f'Splitting {fastafile} in batches of {batch_size} and writing to {prefix}')
    record_iter = SeqIO.parse(open(fastafile, 'r'),"fasta")
    for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
        filename = f"{prefix}_batch_{i + 1}.fasta"
        with open(filename, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, filename))