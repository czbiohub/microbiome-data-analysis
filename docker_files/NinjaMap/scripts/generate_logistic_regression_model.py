#!/usr/bin/env python3

import pandas as pd
import sys

bam_summary_file=sys.argv[1]

# load dataset
bam_summ = pd.read_csv(bam_summary_file, header=True)