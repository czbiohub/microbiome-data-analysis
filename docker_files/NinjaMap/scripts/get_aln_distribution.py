#!/usr/bin/env python3
import gzip
import os
import pysam
import re
import sys

bamfile_name = sys.argv[1]
output_name = sys.argv[2]

truthLabel = os.path.basename(bamfile_name).split('.')[0]

# f = open(output_name, 'w')
f = gzip.open(output_name+".gz", 'wt')
f.write('qname,reference_name,align_len,query_len,aln_cov,quality,perc_id,aln_score,mate_score,mismatches,gap_open,gap_ext,is_dup,is_primary,is_supp,is_truth\n')

bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
for aln in bamfile.fetch(until_eof=True):
    if aln.is_paired and aln.is_proper_pair and not aln.is_secondary:
        ref_name = aln.reference_name.split('_')[0]
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
        is_truth = (ref_name == truthLabel) 
        # perfect alignment over the greatest distance.
        f.write(
            str(aln.qname)+','+
            str(ref_name)+','+
            str(align_len)+','+
            str(query_len)+','+
            str(aln_cov)+','+
            str(quality)+','+
            str(perc_id)+','+
            str(aln_score)+','+
            str(mate_score)+','+
            str(mismatches)+','+
            str(gap_open)+','+
            str(gap_ext)+','+
            str(is_dup)+','+
            str(is_primary)+','+
            str(is_supp)+','+
            str(is_truth)+'\n'
            )