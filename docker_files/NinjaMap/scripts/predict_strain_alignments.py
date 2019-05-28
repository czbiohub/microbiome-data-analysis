#!/usr/bin/env python3
import gzip
import numpy as np
import os
import pysam
import pandas as pd
import re
import sys

bamfile_name = sys.argv[1]
output_name = sys.argv[2]

# truthLabel = os.path.basename(bamfile_name).split('.')[0]

# f = open(output_name, 'w')
# f = gzip.open(output_name+".gz", 'wt')
# f.write('qname,reference_name,align_len,query_len,aln_cov,quality,perc_id,aln_score,mate_score,mismatches,gap_open,gap_ext,is_dup,is_primary,is_supp\n')
df = pd.DataFrame(columns = ['qname','reference_name','align_len','query_len','aln_cov','quality','perc_id','aln_score','mate_score','mismatches','gap_open','gap_ext','is_dup','is_primary','is_supp'])
bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
for aln in bamfile.fetch(until_eof=True):
    if aln.is_paired and aln.is_proper_pair and not aln.is_secondary:
        ref_name = aln.reference_name
        align_len = aln.query_alignment_length
        query_len = aln.query_length
        aln_cov = align_len/float(query_len)
        quality = aln.mapping_quality
        perc_id = 100*(align_len-dict(aln.tags)['NM'])/float(align_len)
        aln_score = dict(aln.tags)['AS']
        mate_score = dict(aln.tags)['YS']
        mismatches = dict(aln.tags)['XM']
        gap_open = dict(aln.tags)['XO']
        gap_ext = dict(aln.tags)['XG']
        is_dup = aln.is_duplicate
        is_primary = not(aln.is_secondary)
        is_supp = aln.is_supplementary
        # is_truth = (ref_name == truthLabel) 
        # perfect alignment over the greatest distance.
        aln_data={
            'qname' : aln.qname,
            'reference_name' : ref_name,
            'align_len' : align_len,
            'query_len' : query_len,
            'aln_cov' : aln_cov,
            'quality' : quality,
            'perc_id' : perc_id,
            'aln_score' : aln_score,
            'mate_score' : mate_score,
            'mismatches' : mismatches,
            'gap_open' : gap_open,
            'gap_ext' : gap_ext,
            'is_dup' : is_dup,
            'is_primary' : is_primary,
            'is_supp' : is_supp
        }
        df.append(aln_data, ignore_index=False)


feature_cols = ["align_len","query_len","aln_cov","quality","perc_id","aln_score","mate_score","mismatches","gap_open","gap_ext","is_dup","is_primary","is_supp"]
# real_test_file = "/Users/sunit.jain/Research/SyntheticCommunities/StrainAbundance/Dorea-longicatena-DSM-13814/Dorea-longicatena-DSM-13814.processed.sortedByCoord.csv"
# real_test = pd.read_csv(real_test_file, header=0)

predict_df = df[feature_cols]
predictions = logreg.predict(predict_df)

predictions.to_csv(output_name,sep=',')

print("Accuracy:",metrics.accuracy_score(real_test_actual, real_test_pred))
print("Precision:",metrics.precision_score(real_test_actual, real_test_pred))
print("Recall:",metrics.recall_score(real_test_actual, real_test_pred))